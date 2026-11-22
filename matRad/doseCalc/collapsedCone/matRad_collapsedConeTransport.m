function dose = matRad_collapsedConeTransport(terma, rho, voxelSize, dirs, weights, att, origin)
% matRad_collapsedConeTransport MATLAB reference implementation of the
% collapsed cone convolution transport (same algorithm as
% matRad_collapsedConeMex)
%
% For each direction the volume is traversed plane by plane along the
% dominant axis on a sheared lattice, so that every voxel is visited exactly
% once per direction. The lattice is anchored to a global grid through the
% origin of the arrays, so that sub-volumes of the same grid give consistent
% results. Along each line the Ahnesjo recursion is applied to every
% exponential kernel component with the radiological step length
% dt = (voxel size / |cos to dominant axis|) * rho.
%
% call:
%   dose = matRad_collapsedConeTransport(terma, rho, voxelSize, dirs, weights, att)
%   dose = matRad_collapsedConeTransport(terma, rho, voxelSize, dirs, weights, att, origin)
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
n = ones(1, 3);
n(1:ndims(terma)) = size(terma);
s = double(voxelSize(:)');
dose = zeros(size(terma));

for d = 1:size(dirs, 2)
    dose = matRad_transportDirection(dose, terma, rho, n, s, double(origin(:)'), double(dirs(:, d)'), ...
                                     double(weights(:, d)), double(att(:, d)));
end

end

function dose = matRad_transportDirection(dose, terma, rho, n, s, origin, direction, weights, att)
[~, k] = max(abs(direction) ./ s);
t1 = mod(k, 3) + 1;
t2 = mod(k + 1, 3) + 1;
stride = [1 n(1) n(1) * n(2)];
nk = n(k);

step = s(k) / abs(direction(k));
sigma1 = direction(t1) * s(k) / (direction(k) * s(t1));
sigma2 = direction(t2) * s(k) / (direction(k) * s(t2));

globalPlane = (0:nk - 1) + origin(k);
sh1 = floor(globalPlane * sigma1 + 0.5) - origin(t1);
sh2 = floor(globalPlane * sigma2 + 0.5) - origin(t2);

[i0, j0] = ndgrid(-max(sh1):n(t1) - 1 - min(sh1), -max(sh2):n(t2) - 1 - min(sh2));
i0 = i0(:);
j0 = j0(:);
flow = zeros(numel(i0), numel(weights));

if direction(k) > 0
    planes = 0:nk - 1;
else
    planes = nk - 1:-1:0;
end

for plane = planes
    i = i0 + sh1(plane + 1);
    j = j0 + sh2(plane + 1);
    valid = find(i >= 0 & i < n(t1) & j >= 0 & j < n(t2));
    if isempty(valid)
        continue
    end
    ix = 1 + plane * stride(k) + i(valid) * stride(t1) + j(valid) * stride(t2);
    dt = step * rho(ix);
    t = terma(ix);
    deposit = zeros(size(ix));
    for c = 1:numel(weights)
        [flow(valid, c), dc] = matRad_depositStep(dt, t, weights(c), att(c), flow(valid, c));
        deposit = deposit + dc;
    end
    dose(ix) = dose(ix) + deposit;
end
end

function [p, d] = matRad_depositStep(dt, t, weight, att, p)
if weight == 0
    d = zeros(size(p));
    return
end
x = att * dt;
oneMinusE = zeros(size(x));
g = ones(size(x));
pos = x > 0;
oneMinusE(pos) = -expm1(-x(pos));
g(pos) = oneMinusE(pos) ./ x(pos);
source = t * weight;
d = p .* g + source .* (1 - g);
p = p .* (1 - oneMinusE) + source .* oneMinusE;
end
