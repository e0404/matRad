function ccKernel = matRad_buildCollapsedConeKernel(runDir, spectrumFile, outFile, numPolarAngles)
% matRad_buildCollapsedConeKernel builds the collapsed cone kernel data file
% from monoenergetic EGSnrc EDKnrc runs
%
% Reads all .egslst files in runDir, combines them into a polyenergetic
% kernel weighted by TERMA, fits the collapsed cone parameters and stores
% the result. The missing energy fraction (escape) is kept and reported,
% no renormalization is applied.
%
% call:
%   ccKernel = matRad_buildCollapsedConeKernel(runDir, spectrumFile)
%   ccKernel = matRad_buildCollapsedConeKernel(runDir, spectrumFile, outFile)
%   ccKernel = matRad_buildCollapsedConeKernel(runDir, spectrumFile, outFile, numPolarAngles)
%
% input:
%   runDir:         folder with the EDKnrc .egslst files
%   spectrumFile:   EGSnrc spectrum file, e.g.
%                   $HEN_HOUSE/spectra/egsnrc/mohan6.spectrum
%   outFile:        (optional) output .mat file. Default:
%                   matRad/doseCalc/collapsedCone/matRad_collapsedConeKernel_6MV.mat
%   numPolarAngles: (optional) number of polar bins of the default fit.
%                   Default: 12
%
% output:
%   ccKernel:       kernel struct as stored in outFile
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

if nargin < 3 || isempty(outFile)
    outFile = fullfile(matRad_cfg.matRadSrcRoot, 'doseCalc', 'collapsedCone', 'matRad_collapsedConeKernel_6MV.mat');
end
if nargin < 4
    numPolarAngles = 12;
end

files = dir(fullfile(runDir, '*.egslst'));
if isempty(files)
    matRad_cfg.dispError('No .egslst files found in %s.', runDir);
end

mono = [];
for i = 1:numel(files)
    mono = [mono matRad_readEdkListing(fullfile(runDir, files(i).name))]; %#ok<AGROW>
end
[~, order] = sort([mono.energy]);
mono = mono(order);

matRad_cfg.dispInfo('Monoenergetic kernels:\n');
matRad_reportEdkStatistics(mono);

[poly, spectrum] = matRad_combineEdkSpectrum(mono, spectrumFile);
matRad_cfg.dispInfo('Polyenergetic kernel (%s):\n', spectrum.title);
polyStats = matRad_reportEdkStatistics(poly);

fit = matRad_fitCollapsedConeKernel(poly, numPolarAngles);
matRad_reportFit(fit);

ccKernel = struct();
ccKernel.meta.description = 'Polyenergetic 6 MV collapsed cone kernel in water';
ccKernel.meta.source = 'EGSnrc EDKnrc, source 0, 521icru.pegs4dat, H2O521ICRU';
ccKernel.meta.spectrum = spectrum.title;
ccKernel.meta.weighting = 'TERMA fraction (fluence * E * mu/rho of water)';
ccKernel.meta.radialUnit = 'mm of water';
ccKernel.meta.numHistoriesPerEnergy = [mono.numHistories];
ccKernel.meta.created = datestr(now, 'yyyy-mm-dd');
ccKernel.missingFraction = polyStats.missing;
ccKernel.cones = rmfield(poly, {'energy', 'numHistories'});
ccKernel.mono = rmfield(mono, 'numHistories');
ccKernel.spectrum = spectrum;
ccKernel.fit = fit;

save(outFile, 'ccKernel', '-v7');
matRad_cfg.dispInfo('Kernel saved to %s\n', outFile);

end

function matRad_reportFit(fit)
matRad_cfg = MatRad_Config.instance();
matRad_cfg.dispInfo('%7s %9s %9s %9s %9s %9s %9s %9s %9s\n', 'theta', 'A', 'a [1/mm]', 'B', 'b [1/mm]', ...
                    'cumErrA', 'cumErrB', 'tailA', 'tailB');
for p = 1:fit.numPolarAngles
    matRad_cfg.dispInfo('%7.2f %9.3e %9.3e %9.3e %9.3e %8.1f%% %8.1f%% %8.2f%% %8.2f%%\n', fit.polarAngles(p), ...
                        fit.A(p), fit.a(p), fit.B(p), fit.b(p), 100 * fit.residual.primaryCum(p), ...
                        100 * fit.residual.scatterCum(p), 100 * fit.residual.primaryEnergy(p), ...
                        100 * fit.residual.scatterEnergy(p));
end
matRad_cfg.dispInfo('sum(A) = %.4f, sum(B) = %.4f, sum(A+B) = %.4f, missing = %.4f\n', sum(fit.A), sum(fit.B), ...
                    sum(fit.A + fit.B), 1 - sum(fit.A + fit.B));
end
