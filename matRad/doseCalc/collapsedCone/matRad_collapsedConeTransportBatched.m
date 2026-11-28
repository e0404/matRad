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
% The transport is computed in the precision of the TERMA: a single TERMA
% gives a single dose (deviation from the double result below 1e-6 of the
% maximum dose), any other type is computed in double. If the TERMA or the
% density is a gpuArray, the transport runs on the GPU and the dose is
% returned as gpuArray. The lattice is always computed in double on the
% host, so that it is the same for all precisions.
%
% call:
%   dose = matRad_collapsedConeTransportBatched(terma, rho, voxelSize, dirs, weights, att)
%   dose = matRad_collapsedConeTransportBatched(terma, rho, voxelSize, dirs, weights, att, origin)
%
% input:
%   terma:      3D array of TERMA (double or single, optionally gpuArray)
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
%   dose:       3D array of deposited dose (same size, precision and
%               device as terma)
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

% working precision and device follow the TERMA (and the density for the device)
precision = 'double';
if strcmp(matRad_underlyingTypeCompat(terma), 'single')
    precision = 'single';
end
useGPU = matRad_isgpuarrayCompat(terma) || matRad_isgpuarrayCompat(rho);
terma = matRad_toWorkArray(terma, precision, useGPU);
rho = matRad_toWorkArray(rho, precision, useGPU);

% geometry and kernel parameters stay in double on the host
dirs = double(matRad_gatherCompat(dirs));
weights = double(matRad_gatherCompat(weights));
att = double(matRad_gatherCompat(att));
origin = double(matRad_gatherCompat(origin(:)'));
n = ones(1, 3);
n(1:ndims(terma)) = size(terma);
s = double(matRad_gatherCompat(voxelSize(:)'));
numDirs = size(dirs, 2);
dose = matRad_zerosLike(n, terma);

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
    dosePerm = matRad_zerosLike(geometry.n, terma);
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
precision = matRad_underlyingTypeCompat(terma);
useGPU = matRad_isgpuarrayCompat(terma);
weights = matRad_toWorkArray(reshape(weights(keep), 1, 1, []), precision, useGPU);
negAttStep = matRad_toWorkArray(reshape(-attStep(keep), 1, 1, []), precision, useGPU);
numColumns = numel(dirOfColumn);
sumWeights = sum(weights);

% rows and columns of each plane that contain TERMA and planes with voxels
% without density (on the host, so that the loop needs no device queries)
termaRows = matRad_gatherCompat(reshape(any(terma, 2), n1, nk));
termaCols = matRad_gatherCompat(reshape(any(terma, 1), n2, nk));
hasVoid = matRad_gatherCompat(reshape(any(any(~(rho > 0), 1), 2), nk, 1));

% the flow carries a border of one voxel that stays zero, so that shifting
% moves zeros in at the edges of the volume
flow = matRad_zerosLike([n1 + 2, n2 + 2, numColumns], terma);
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

    % Ahnesjo recursion for all columns, see matRad_collapsedConeTransport,
    % with negExpm1 = -(1 - exp(-x)) and g = (1 - exp(-x)) / x
    rhoPlane = rho(i1, i2, plane);
    termaPlane = terma(i1, i2, plane);
    flowBox = flow(i1 + 1, i2 + 1, :);
    if useGPU
        % one fused element-wise kernel instead of many small array operations
        [deposit, flowBox] = arrayfun(@matRad_recursionStep, flowBox, termaPlane, rhoPlane, negAttStep, weights);
    else
        negX = rhoPlane .* negAttStep;
        negExpm1 = expm1(negX);
        g = negExpm1 ./ negX;
        if hasVoid(plane)
            g(repmat(~(rhoPlane > 0), [1 1 numColumns])) = 1;
        end
        flowMinusSource = flowBox - termaPlane .* weights;
        deposit = flowMinusSource .* g;
        flowBox = flowBox + flowMinusSource .* negExpm1;
    end
    % deposit = flow .* g + source .* (1 - g), new flow = flow .* exp(-x) + source .* (1 - exp(-x))
    dose(i1, i2, plane) = dose(i1, i2, plane) + termaPlane .* sumWeights + sum(deposit, 3);
    flow(i1 + 1, i2 + 1, :) = flowBox;
end
end

function [deposit, flow] = matRad_recursionStep(flow, terma, rho, negAttStep, weight)
% Scalar Ahnesjo recursion step for gpuArray arrayfun; has to stay in line
% with the array expressions in matRad_transportGroup. The deposit lacks
% the term terma * weight, which the caller adds for all columns at once.
negX = rho * negAttStep;
negExpm1 = expm1(negX);
if rho > 0
    g = negExpm1 / negX;
else
    g = negExpm1 - negExpm1 + 1; % one in the working precision
end
flowMinusSource = flow - terma * weight;
deposit = flowMinusSource * g;
flow = flow + flowMinusSource * negExpm1;
end

function [flow, box] = matRad_shiftFlow(flow, d1, d2, box)
% Shifts the flow of every column by d1 and d2 voxels (-1, 0 or 1) along the
% two transverse axes. box is the bounding box of the non-zero flow in voxel
% indices; the flow array has a zero border of one voxel (index offset 1).
% Flow shifted out of the volume is discarded. Zeros are shifted in, since
% lattice lines enter the volume without flow: the destination covers the
% box plus one voxel, and its source is zero outside the box.
n1 = size(flow, 1) - 2;
n2 = size(flow, 2) - 2;
rows = max(box(1) - 1, 1):min(box(2) + 1, n1);
cols = max(box(3) - 1, 1):min(box(4) + 1, n2);

shiftCode = 3 * d1 + d2;
codes = unique(shiftCode(shiftCode ~= 0));
for q = 1:numel(codes)
    columns = find(shiftCode == codes(q));
    a = d1(columns(1));
    b = d2(columns(1));
    flow(rows + 1, cols + 1, columns) = flow(rows + 1 - a, cols + 1 - b, columns);
end

box = [max(box(1) + min(d1), 1), min(box(2) + max(d1), n1), max(box(3) + min(d2), 1), min(box(4) + max(d2), n2)];
end

function x = matRad_toWorkArray(x, precision, useGPU)
% Casts x to the working precision and moves it to the GPU if requested
x = cast(x, precision);
if useGPU && ~matRad_isgpuarrayCompat(x)
    x = gpuArray(x);
end
end

function z = matRad_zerosLike(sz, prototype)
% Zeros of the precision and on the device of prototype
if matRad_isgpuarrayCompat(prototype)
    z = zeros(sz, 'like', prototype);
else
    z = zeros(sz, class(prototype));
end
end
