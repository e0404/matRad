function kernel = matRad_fitCollapsedConeKernel(edk, numPolarAngles, minRadius, numScatterExponentials)
% matRad_fitCollapsedConeKernel collapses a conical energy deposition kernel
% into polar bins and fits a sum of exponentials per component
%
% For each polar bin the energy fraction deposited per unit radial distance
% r (mm of water) is modeled as
%   dE/dr = A a exp(-a r) + sum_k B_k b_k exp(-b_k r)
% where (A,a) describe the primary and (B_k,b_k) the scatter component with
% numScatterExponentials terms. A and sum_k B_k are the primary and scatter
% energy fractions of the bin found in the Monte Carlo data (no
% renormalization, escaping energy stays missing). The attenuation
% coefficients and the split of the scatter energy between its terms are
% fitted to the normalized cumulative radial energy distribution of each
% component.
%
% call:
%   kernel = matRad_fitCollapsedConeKernel(edk, numPolarAngles)
%   kernel = matRad_fitCollapsedConeKernel(edk, numPolarAngles, minRadius)
%   kernel = matRad_fitCollapsedConeKernel(edk, numPolarAngles, minRadius, numScatterExponentials)
%
% input:
%   edk:                    conical kernel struct with fields angleEdges
%                           [deg], radiusEdges [mm], primary and scatter
%                           (nCones x nShells)
%   numPolarAngles:         number of polar bins, must divide the number of
%                           cones
%   minRadius:              (optional) smallest radius [mm] included in the
%                           fit and in the residuals. Default: 0
%   numScatterExponentials: (optional) number of exponential terms of the
%                           scatter component (1 or 2). Default: 1
%
% output:
%   kernel:         struct with fields
%                     numPolarAngles, minRadius, numScatterExponentials
%                     polarEdges      bin edges [deg]
%                     polarAngles     solid angle bisector of each bin [deg]
%                     A, a            primary weight and attenuation
%                                     coefficient (1 x nBins, 1/mm of water)
%                     B, b            scatter weights and attenuation
%                                     coefficients (nTerms x nBins)
%                     residual        struct with, per bin, the maximum
%                                     absolute error of the normalized
%                                     cumulative energy over all radii
%                                     (primaryCum, scatterCum), for
%                                     r <= 10 mm (primaryCumNear,
%                                     scatterCumNear) and for r >= 100 mm
%                                     (primaryCumFar, scatterCumFar), and the
%                                     fraction of the modeled energy beyond
%                                     the kernel range (primaryEnergy,
%                                     scatterEnergy, negative)
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

if nargin < 3 || isempty(minRadius)
    minRadius = 0;
end
if nargin < 4 || isempty(numScatterExponentials)
    numScatterExponentials = 1;
end
if ~any(numScatterExponentials == [1 2])
    matRad_cfg.dispError('Only one or two scatter exponentials are supported.');
end

nCones = size(edk.primary, 1);
if mod(nCones, numPolarAngles) ~= 0
    matRad_cfg.dispError('Number of polar angles (%d) must divide the number of kernel cones (%d).', numPolarAngles, nCones);
end
conesPerBin = nCones / numPolarAngles;

kernel = struct();
kernel.numPolarAngles = numPolarAngles;
kernel.minRadius = minRadius;
kernel.numScatterExponentials = numScatterExponentials;
kernel.polarEdges = edk.angleEdges(1:conesPerBin:end);
cosEdges = cosd(kernel.polarEdges);
kernel.polarAngles = acosd(0.5 * (cosEdges(1:end - 1) + cosEdges(2:end)));

rOuter = edk.radiusEdges(2:end);
comps = {'primary', 'scatter'};
numTerms = [1 numScatterExponentials];
inFit = rOuter >= minRadius;
masks = {true(size(rOuter)), inFit & rOuter <= 10, inFit & rOuter >= 100};
suffix = {'Cum', 'CumNear', 'CumFar'};

for c = 1:2
    val = edk.(comps{c});
    weight = zeros(numTerms(c), numPolarAngles);
    att = ones(numTerms(c), numPolarAngles);
    err = zeros(3, numPolarAngles);
    tail = zeros(1, numPolarAngles);
    for p = 1:numPolarAngles
        rows = (p - 1) * conesPerBin + (1:conesPerBin);
        y = sum(val(rows, :), 1);
        energy = sum(y);
        if energy <= 0
            continue
        end
        cumData = cumsum(y) / energy;
        [frac, att(:, p)] = matRad_fitCumulativeExponentials(rOuter(inFit), cumData(inFit), numTerms(c));
        weight(:, p) = energy * frac;
        model = 1 - frac' * exp(-att(:, p) * rOuter);
        for m = 1:3
            if any(masks{m})
                err(m, p) = max(abs(model(masks{m}) - cumData(masks{m})));
            end
        end
        tail(p) = -frac' * exp(-att(:, p) * rOuter(end));
    end
    if c == 1
        kernel.A = weight;
        kernel.a = att;
    else
        kernel.B = weight;
        kernel.b = att;
    end
    for m = 1:3
        kernel.residual.([comps{c} suffix{m}]) = err(m, :);
    end
    kernel.residual.([comps{c} 'Energy']) = tail;
end

end

function [frac, att] = matRad_fitCumulativeExponentials(r, cumData, numTerms)
r = r(:)';
cumData = cumData(:)';
if numTerms == 1
    cost = @(logAtt) sum((1 - exp(-exp(logAtt) * r) - cumData).^2);
    logGrid = linspace(log(1e-5), log(20), 300);
    [~, iBest] = min(arrayfun(cost, logGrid));
    lo = logGrid(max(iBest - 1, 1));
    hi = logGrid(min(iBest + 1, numel(logGrid)));
    att = exp(fminbnd(cost, lo, hi, optimset('TolX', 1e-10, 'Display', 'off')));
    frac = 1;
    return
end

logGrid = linspace(log(1e-5), log(20), 80);
bestCost = Inf;
bestLog = [0 0];
for i = 1:numel(logGrid)
    for j = 1:i - 1
        c = matRad_twoTermCost(logGrid([i j]), r, cumData);
        if c < bestCost
            bestCost = c;
            bestLog = logGrid([i j]);
        end
    end
end
logAtt = fminsearch(@(la) matRad_twoTermCost(la, r, cumData), bestLog, ...
                    optimset('TolX', 1e-10, 'TolFun', 1e-14, 'MaxFunEvals', 4000, 'MaxIter', 4000, 'Display', 'off'));
[~, f1] = matRad_twoTermCost(logAtt, r, cumData);
att = exp(logAtt(:));
frac = [f1; 1 - f1];
[att, order] = sort(att, 'descend');
frac = frac(order);
end

function [cost, f1] = matRad_twoTermCost(logAtt, r, cumData)
e1 = exp(-exp(logAtt(1)) * r);
e2 = exp(-exp(logAtt(2)) * r);
target = 1 - cumData - e2;
basis = e1 - e2;
f1 = sum(target .* basis) / max(sum(basis.^2), realmin);
f1 = min(max(f1, 0), 1);
cost = sum((f1 * e1 + (1 - f1) * e2 - (1 - cumData)).^2);
end
