function kernel = matRad_fitCollapsedConeKernel(edk, numPolarAngles, minRadius)
% matRad_fitCollapsedConeKernel collapses a conical energy deposition kernel
% into polar bins and fits the Ahnesjo two-exponential parametrization
%
% For each polar bin the energy fraction deposited per unit radial distance
% r (mm of water) is modeled as
%   dE/dr = A a exp(-a r) + B b exp(-b r)
% where (A,a) describe the primary and (B,b) the scatter component.
% A and B are the energy fractions of the bin found in the Monte Carlo data
% (no renormalization, escaping energy stays missing). The attenuation
% coefficients a and b are fitted to the normalized cumulative radial
% energy distribution of each component.
%
% call:
%   kernel = matRad_fitCollapsedConeKernel(edk, numPolarAngles)
%   kernel = matRad_fitCollapsedConeKernel(edk, numPolarAngles, minRadius)
%
% input:
%   edk:            conical kernel struct with fields angleEdges [deg],
%                   radiusEdges [mm], primary and scatter (nCones x nShells)
%   numPolarAngles: number of polar bins, must divide the number of cones
%   minRadius:      (optional) smallest radius [mm] included in the fit
%                   and in the residuals. Default: 0
%
% output:
%   kernel:         struct with fields
%                     numPolarAngles
%                     polarEdges      bin edges [deg]
%                     polarAngles     solid angle bisector of each bin [deg]
%                     A, a, B, b      fit parameters (a, b in 1/mm of water)
%                     residual        struct with, per bin, the maximum
%                                     absolute error of the normalized
%                                     cumulative energy over all radii
%                                     (primaryCum, scatterCum) and for
%                                     r <= 10 mm (primaryCumNear,
%                                     scatterCumNear), and the relative
%                                     error of the energy within the kernel
%                                     range (primaryEnergy, scatterEnergy)
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

matRad_cfg = MatRad_Config.instance();

if nargin < 3
    minRadius = 0;
end

nCones = size(edk.primary, 1);
if mod(nCones, numPolarAngles) ~= 0
    matRad_cfg.dispError('Number of polar angles (%d) must divide the number of kernel cones (%d).', numPolarAngles, nCones);
end
conesPerBin = nCones / numPolarAngles;

kernel = struct();
kernel.numPolarAngles = numPolarAngles;
kernel.minRadius = minRadius;
kernel.polarEdges = edk.angleEdges(1:conesPerBin:end);
cosEdges = cosd(kernel.polarEdges);
kernel.polarAngles = acosd(0.5 * (cosEdges(1:end - 1) + cosEdges(2:end)));

rOuter = edk.radiusEdges(2:end);
comps = {'primary', 'scatter'};
amp = zeros(2, numPolarAngles);
att = ones(2, numPolarAngles);
cumErr = zeros(2, numPolarAngles);
cumErrNear = zeros(2, numPolarAngles);
inFit = rOuter >= minRadius;
isNear = inFit & rOuter <= 10;
energyErr = zeros(2, numPolarAngles);

for c = 1:2
    val = edk.(comps{c});
    for p = 1:numPolarAngles
        rows = (p - 1) * conesPerBin + (1:conesPerBin);
        y = sum(val(rows, :), 1);
        amp(c, p) = sum(y);
        if amp(c, p) <= 0
            continue
        end
        cumData = cumsum(y) / amp(c, p);
        [att(c, p), cumErr(c, p)] = matRad_fitCumulativeExponential(rOuter(inFit), cumData(inFit));
        cumErrNear(c, p) = max(abs(1 - exp(-att(c, p) * rOuter(isNear)) - cumData(isNear)));
        energyErr(c, p) = -exp(-att(c, p) * rOuter(end));
    end
end

kernel.A = amp(1, :);
kernel.a = att(1, :);
kernel.B = amp(2, :);
kernel.b = att(2, :);
kernel.residual.primaryCum = cumErr(1, :);
kernel.residual.scatterCum = cumErr(2, :);
kernel.residual.primaryCumNear = cumErrNear(1, :);
kernel.residual.scatterCumNear = cumErrNear(2, :);
kernel.residual.primaryEnergy = energyErr(1, :);
kernel.residual.scatterEnergy = energyErr(2, :);

end

function [att, maxErr] = matRad_fitCumulativeExponential(r, cumData)
cost = @(logAtt) sum((1 - exp(-exp(logAtt) * r) - cumData).^2);
logGrid = linspace(log(1e-5), log(20), 300);
costGrid = arrayfun(cost, logGrid);
[~, iBest] = min(costGrid);
lo = logGrid(max(iBest - 1, 1));
hi = logGrid(min(iBest + 1, numel(logGrid)));
att = exp(fminbnd(cost, lo, hi, optimset('TolX', 1e-10, 'Display', 'off')));
maxErr = max(abs(1 - exp(-att * r) - cumData));
end
