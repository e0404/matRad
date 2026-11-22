function stats = matRad_reportEdkStatistics(edk, numOuterShells)
% matRad_reportEdkStatistics reports energy fractions and statistical
% uncertainties of conical energy deposition kernels
%
% call:
%   stats = matRad_reportEdkStatistics(edk)
%   stats = matRad_reportEdkStatistics(edk, numOuterShells)
%
% input:
%   edk:            struct array of kernels as returned by
%                   matRad_readEdkListing or matRad_combineEdkSpectrum
%   numOuterShells: (optional) number of outermost radial shells for the
%                   scatter uncertainty report. Default: 3
%
% output:
%   stats:          struct array with energy, total, primary, scatter and
%                   missing energy fractions, forward/backward ratios and
%                   relative scatter uncertainties of the outer shells
%                   (angle integrated and worst cone)
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

if nargin < 2
    numOuterShells = 3;
end

stats = struct([]);
for i = 1:numel(edk)
    k = edk(i);
    nShells = size(k.total, 2);
    outer = nShells - numOuterShells + 1:nShells;
    forward = k.angleEdges(2:end) <= 90;

    s.energy = k.energy;
    s.total = sum(k.total(:));
    s.primary = sum(k.primary(:));
    s.scatter = sum(k.scatter(:));
    s.missing = 1 - s.total;
    s.primaryForwardBackward = sum(sum(k.primary(forward, :))) / sum(sum(k.primary(~forward, :)));
    s.scatterForwardBackward = sum(sum(k.scatter(forward, :))) / sum(sum(k.scatter(~forward, :)));
    s.outerRadii = k.radiusEdges(outer + 1);

    absUnc = k.scatterRelUnc .* k.scatter;
    absUnc(~isfinite(absUnc)) = 0;
    shellSum = sum(k.scatter(:, outer), 1);
    s.outerScatterRelUnc = sqrt(sum(absUnc(:, outer).^2, 1)) ./ shellSum;
    s.outerScatterMaxConeRelUnc = max(k.scatterRelUnc(:, outer), [], 1);
    stats = [stats s]; %#ok<AGROW>
end

matRad_cfg.dispInfo('%8s %8s %8s %8s %8s %9s %9s', 'E [MeV]', 'total', 'primary', 'scatter', 'missing', 'pri f/b', 'sca f/b');
for j = 1:numOuterShells
    matRad_cfg.dispInfo(' sca@%3.0fcm', stats(1).outerRadii(j) / 10);
end
matRad_cfg.dispInfo(' (angle-integrated / worst cone)\n');
for i = 1:numel(stats)
    s = stats(i);
    matRad_cfg.dispInfo('%8.3f %8.4f %8.4f %8.4f %8.4f %9.1f %9.2f', s.energy, s.total, s.primary, s.scatter, s.missing, ...
                        s.primaryForwardBackward, s.scatterForwardBackward);
    for j = 1:numOuterShells
        matRad_cfg.dispInfo(' %4.1f/%4.0f%%', 100 * s.outerScatterRelUnc(j), 100 * s.outerScatterMaxConeRelUnc(j));
    end
    matRad_cfg.dispInfo('\n');
end

end
