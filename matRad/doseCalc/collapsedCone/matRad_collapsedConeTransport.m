function dose = matRad_collapsedConeTransport(terma, rho, voxelSize, dirs, wA, aA, wB, aB)
% matRad_collapsedConeTransport MATLAB reference implementation of the
% collapsed cone convolution transport (same algorithm as
% matRad_collapsedConeMex)
%
% For each direction the volume is traversed plane by plane along the
% dominant axis on a sheared lattice, so that every voxel is visited exactly
% once per direction. Along each line the Ahnesjo recursion is applied to
% the primary (A,a) and scatter (B,b) components with the radiological step
% length dt = (voxel size / |cos to dominant axis|) * rho.
%
% call:
%   dose = matRad_collapsedConeTransport(terma, rho, voxelSize, dirs, A, a, B, b)
%
% input:
%   terma:      3D array of TERMA
%   rho:        3D array of relative electron densities (same size)
%   voxelSize:  voxel size along the three array dimensions [mm]
%   dirs:       3 x nDirs unit transport directions in array axes
%   A, a, B, b: 1 x nDirs kernel weights and attenuation coefficients
%               [1/mm of water] per direction
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

terma = double(terma);
rho = double(rho);
n = ones(1, 3);
n(1:ndims(terma)) = size(terma);
s = double(voxelSize(:)');
dose = zeros(size(terma));

for d = 1:size(dirs, 2)
    dose = matRad_transportDirection(dose, terma, rho, n, s, double(dirs(:, d)'), ...
                                     [wA(d) aA(d)], [wB(d) aB(d)]);
end

end

function dose = matRad_transportDirection(dose, terma, rho, n, s, direction, ca, cb)
[~, k] = max(abs(direction) ./ s);
t1 = mod(k, 3) + 1;
t2 = mod(k + 1, 3) + 1;
stride = [1 n(1) n(1) * n(2)];

step = s(k) / abs(direction(k));
delta1 = direction(t1) * s(k) / (abs(direction(k)) * s(t1));
delta2 = direction(t2) * s(k) / (abs(direction(k)) * s(t2));
nk = n(k);

q = 0:nk - 1;
sh1 = floor(q * delta1 + 0.5);
sh2 = floor(q * delta2 + 0.5);

[i0, j0] = ndgrid(-max(sh1(end), 0):n(t1) - 1 - min(sh1(end), 0), ...
                  -max(sh2(end), 0):n(t2) - 1 - min(sh2(end), 0));
i0 = i0(:);
j0 = j0(:);
pa = zeros(size(i0));
pb = zeros(size(i0));

for iq = 1:nk
    i = i0 + sh1(iq);
    j = j0 + sh2(iq);
    valid = find(i >= 0 & i < n(t1) & j >= 0 & j < n(t2));
    if isempty(valid)
        continue
    end
    if direction(k) > 0
        plane = iq - 1;
    else
        plane = nk - iq;
    end
    ix = 1 + plane * stride(k) + i(valid) * stride(t1) + j(valid) * stride(t2);
    dt = step * rho(ix);
    t = terma(ix);
    [pa(valid), da] = matRad_depositStep(dt, t, ca, pa(valid));
    [pb(valid), db] = matRad_depositStep(dt, t, cb, pb(valid));
    dose(ix) = dose(ix) + (da + db);
end
end

function [p, d] = matRad_depositStep(dt, t, c, p)
if c(1) == 0
    d = zeros(size(p));
    return
end
x = c(2) * dt;
oneMinusE = zeros(size(x));
g = ones(size(x));
pos = x > 0;
oneMinusE(pos) = -expm1(-x(pos));
g(pos) = oneMinusE(pos) ./ x(pos);
source = t * c(1);
d = p .* g + source .* (1 - g);
p = p .* (1 - oneMinusE) + source .* oneMinusE;
end
