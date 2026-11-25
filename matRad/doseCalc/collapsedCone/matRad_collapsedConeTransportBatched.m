function dose = matRad_collapsedConeTransportBatched(terma, rho, voxelSize, dirs, weights, att, origin)
% matRad_collapsedConeTransportBatched vectorized MATLAB implementation of
% the collapsed cone convolution transport (same algorithm and lattice as
% matRad_collapsedConeTransport and matRad_collapsedConeMex)
%
% matRad_collapsedConeTransport loops over the directions and gathers the
% voxels of the lattice lines plane by plane. Here all directions (and
% kernel components) that share the dominant axis and its sign are
% transported together: their flow is stored aligned to the voxels of the
% current plane as an n1 x n2 x nColumns array (one column per direction and
% component), which is shifted by at most one voxel per transverse axis
% whenever the lattice shift of a direction changes between two planes.
% A plane is then processed with slicing, implicit expansion and a sum over
% the columns only, i.e. without index gather/scatter, in six plane sweeps
% instead of one per direction. Planes and voxels the flow has not reached
% yet are skipped by tracking the bounding box of the non-zero flow.
%
% The result equals the one of matRad_collapsedConeTransport up to the
% floating point summation order.
%
% call:
%   dose = matRad_collapsedConeTransportBatched(terma, rho, voxelSize, dirs, weights, att)
%   dose = matRad_collapsedConeTransportBatched(terma, rho, voxelSize, dirs, weights, att, origin)
%
% input:
%   terma:      3D array of TERMA
%   rho:        3D array of relative electron densities (same size)
%   voxelSize:  voxel size along the three array dimensions [mm]
%   dirs:       3 x nDirs unit transport directions in array axes
%   weights:    nComponents x nDirs kernel weights of the exponential
%               components (e.g. primary and scatter) per direction
%   att:        nComponents x nDirs attenuation coefficients [1/mm of water]
%   origin:     (optional) zero based index of the first voxel of the arrays
%               in the global grid. Default: [0 0 0]
%
% output:
%   dose:       3D array of deposited dose (same size as terma)
%
% References
%   [1] Ahnesjo, Med. Phys. 16, 577 (1989)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 7 || isempty(origin)
    origin = [0 0 0];
end

matRad_cfg = MatRad_Config.instance();
if size(weights, 2) ~= size(dirs, 2) || ~isequal(size(weights), size(att))
    matRad_cfg.dispError('Kernel weights and attenuation coefficients must be nComponents x nDirs arrays.');
end

terma = double(terma);
rho = double(rho);
dirs = double(dirs);
weights = double(weights);
att = double(att);
origin = double(origin(:)');
n = ones(1, 3);
n(1:ndims(terma)) = size(terma);
s = double(voxelSize(:)');
numDirs = size(dirs, 2);
dose = zeros(n);

% dominant axis of each direction and the sign along it
[~, dominantAxis] = max(abs(dirs) ./ repmat(s(:), 1, numDirs), [], 1);
isForward = dirs(sub2ind(size(dirs), dominantAxis, 1:numDirs)) > 0;

for k = 1:3
    if ~any(dominantAxis == k)
        continue
    end
    % transverse axes first, dominant axis last, so that planes are slices
    t1 = mod(k, 3) + 1;
    t2 = mod(k + 1, 3) + 1;
    perm = [t1 t2 k];
    geometry = struct('n', n(perm), 's', s(perm), 'origin', origin(perm));
    termaPerm = permute(terma, perm);
    rhoPerm = permute(rho, perm);
    dosePerm = zeros(geometry.n);
    for forward = [true false]
        sel = dominantAxis == k & isForward == forward;
        if any(sel)
            dosePerm = matRad_transportGroup(dosePerm, termaPerm, rhoPerm, geometry, dirs(perm, sel), ...
                                             weights(:, sel), att(:, sel), forward);
        end
    end
    dose = dose + ipermute(dosePerm, perm);
end

end

function dose = matRad_transportGroup(dose, terma, rho, geometry, dirs, weights, att, forward)
% Transports all directions of a group. The arrays are permuted such that
% the dominant axis of the group is the third dimension; dirs is permuted
% accordingly (rows: transverse axis 1, transverse axis 2, dominant axis).
n1 = geometry.n(1);
n2 = geometry.n(2);
nk = geometry.n(3);
s = geometry.s;
numComps = size(weights, 1);
numDirs = size(dirs, 2);

% lattice shifts per plane and direction (as in matRad_collapsedConeTransport;
% the constant offset by the origin of the transverse axes does not matter,
% because only the change of the shifts between planes is used)
step = s(3) ./ abs(dirs(3, :));
sigma1 = dirs(1, :) * s(3) ./ (dirs(3, :) * s(1));
sigma2 = dirs(2, :) * s(3) ./ (dirs(3, :) * s(2));
globalPlane = (0:nk - 1)' + geometry.origin(3);
sh1 = floor(globalPlane * sigma1 + 0.5);
sh2 = floor(globalPlane * sigma2 + 0.5);

% one column per direction and component with non-zero weight
attStep = att .* repmat(step, numComps, 1);
dirOfColumn = repmat(1:numDirs, numComps, 1);
keep = weights ~= 0;
if ~any(keep(:))
    return
end
dirOfColumn = dirOfColumn(keep)';
weights = reshape(weights(keep), 1, 1, []);
attStep = reshape(attStep(keep), 1, 1, []);
numColumns = numel(weights);
sumWeights = sum(weights);

% rows and columns of each plane that contain TERMA
termaRows = reshape(any(terma, 2), n1, nk);
termaCols = reshape(any(terma, 1), n2, nk);

flow = zeros(n1, n2, numColumns);
box = [n1 + 1, 0, n2 + 1, 0]; % bounding box [lo1 hi1 lo2 hi2] of the non-zero flow, empty at the start

if forward
    planes = 1:nk;
else
    planes = nk:-1:1;
end

prevPlane = 0;
for plane = planes
    % move the flow to the voxels the lattice lines hit in this plane
    if prevPlane > 0 && box(2) >= box(1)
        d1 = sh1(plane, :) - sh1(prevPlane, :);
        d2 = sh2(plane, :) - sh2(prevPlane, :);
        if any(d1) || any(d2)
            [flow, box] = matRad_shiftFlow(flow, d1(dirOfColumn), d2(dirOfColumn), box);
        end
    end
    prevPlane = plane;

    rows = find(termaRows(:, plane));
    if ~isempty(rows)
        cols = find(termaCols(:, plane));
        box = [min(box(1), rows(1)), max(box(2), rows(end)), min(box(3), cols(1)), max(box(4), cols(end))];
    end
    if box(2) < box(1)
        continue
    end
    i1 = box(1):box(2);
    i2 = box(3):box(4);

    % Ahnesjo recursion for all columns, see matRad_collapsedConeTransport
    rhoPlane = rho(i1, i2, plane);
    termaPlane = terma(i1, i2, plane);
    x = rhoPlane .* attStep;
    oneMinusE = -expm1(-x);
    g = oneMinusE ./ x;
    if ~all(rhoPlane(:) > 0)
        g(repmat(~(rhoPlane > 0), [1 1 numColumns])) = 1;
    end
    flowMinusSource = flow(i1, i2, :) - termaPlane .* weights;
    % deposit = flow .* g + source .* (1 - g), new flow = flow .* (1 - oneMinusE) + source .* oneMinusE
    dose(i1, i2, plane) = dose(i1, i2, plane) + termaPlane .* sumWeights + sum(flowMinusSource .* g, 3);
    flow(i1, i2, :) = flow(i1, i2, :) - flowMinusSource .* oneMinusE;
end
end

function [flow, box] = matRad_shiftFlow(flow, d1, d2, box)
% Shifts the flow of every column by d1 and d2 voxels (-1, 0 or 1) along the
% two transverse axes. Flow shifted out of the array is discarded and zeros
% are shifted in, since lattice lines enter the volume without flow.
n1 = size(flow, 1);
n2 = size(flow, 2);
lo1 = box(1);
hi1 = box(2);
lo2 = box(3);
hi2 = box(4);

shiftCode = 3 * d1 + d2;
codes = unique(shiftCode(shiftCode ~= 0));
for q = 1:numel(codes)
    columns = find(shiftCode == codes(q));
    a = d1(columns(1));
    b = d2(columns(1));
    sourceRows = max(lo1, 1 - a):min(hi1, n1 - a);
    sourceCols = max(lo2, 1 - b):min(hi2, n2 - b);
    block = flow(sourceRows, sourceCols, columns);
    flow(lo1:hi1, lo2:hi2, columns) = 0;
    flow(sourceRows + a, sourceCols + b, columns) = block;
end

box = [max(lo1 + min(d1), 1), min(hi1 + max(d1), n1), max(lo2 + min(d2), 1), min(hi2 + max(d2), n2)];
end
