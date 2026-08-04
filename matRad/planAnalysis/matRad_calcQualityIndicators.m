function qi = matRad_calcQualityIndicators(cst, pln, doseCube, refGy, refVol)
% matRad_calcQualityIndicators calculates DVH quality indicators.
%
% call
%   qi = matRad_calcQualityIndicators(cst,pln,doseCube)
%   qi = matRad_calcQualityIndicators(cst,pln,doseCube,refGy,refVol)
%
% input
%   cst:                matRad cst cell array
%   pln:                matRad pln struct
%   doseCube:           per-fraction dose cube
%   refGy: (optional)   per-fraction dose values for V_XGy calculation
%   refVol:(optional)   volume percentages for D_X calculation
%
% output
%   qi:                 quality indicators in per-fraction dose units
%
% References
%   van't Riet et. al., IJROBP, 1997 Feb 1;37(3):731-6.
%   Kataria et. al., J Med Phys. 2012 Oct-Dec; 37(4)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016-2026 the matRad development team.
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

if ~exist('refGy', 'var')
    refGy = [];
end

if ~exist('refVol', 'var') || isempty(refVol)
    refVol = [2 5 50 95 98];
end

if isempty(refGy)
    finiteDose = doseCube(isfinite(doseCube));
    if isempty(finiteDose)
        refGy = 0;
    else
        refGy = floor(linspace(0, max(finiteDose(:)), 6) * 10) / 10;
    end
end

qi = struct;
targetDoseInfo = matRad_getTargetReferenceDoses(cst, pln);
targetDoseCstIndex = [targetDoseInfo.cstIndex];

for runVoi = 1:size(cst, 1)
    indices = cst{runVoi, 4}{1};
    numOfVoxels = numel(indices);
    voiPrint = sprintf('%3d %20s', cst{runVoi, 1}, cst{runVoi, 2});
    qi(runVoi).name = cst{runVoi, 2};

    doseInVoi = sort(doseCube(indices));

    if ~isempty(doseInVoi)
        qi(runVoi).mean = mean(doseInVoi);
        qi(runVoi).std  = std(doseInVoi);
        qi(runVoi).max  = doseInVoi(end);
        qi(runVoi).min  = doseInVoi(1);

        voiPrint = sprintf(['%s - Mean dose = %5.2f Gy +/- %5.2f Gy ', ...
                            '(Max dose = %5.2f Gy, Min dose = %5.2f Gy)\n%27s'], ...
                           voiPrint, qi(runVoi).mean, qi(runVoi).std, ...
                           qi(runVoi).max, qi(runVoi).min, ' ');

        DX = @(x) matRad_interp1(linspace(0, 1, numOfVoxels), doseInVoi, (100 - x) * 0.01);
        VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;

        for runDX = 1:numel(refVol)
            qi(runVoi).(strcat('D_', num2str(refVol(runDX)))) = DX(refVol(runDX));
            voiPrint = sprintf('%sD%d%% = %5.2f Gy, ', voiPrint, ...
                               refVol(runDX), DX(refVol(runDX)));
        end
        voiPrint = sprintf('%s\n%27s', voiPrint, ' ');

        for runVX = 1:numel(refGy)
            sRefGy = num2str(refGy(runVX), 3);
            qi(runVoi).(['V_' strrep(sRefGy, '.', '_') 'Gy']) = VX(refGy(runVX));
            voiPrint = sprintf(['%sV' sRefGy 'Gy = %6.2f%%, '], ...
                               voiPrint, VX(refGy(runVX)) * 100);
        end
        voiPrint = sprintf('%s\n%27s', voiPrint, ' ');

        if strcmp(cst{runVoi, 3}, 'TARGET') > 0
            targetDoseIx = find(targetDoseCstIndex == runVoi, 1, 'first');
            if isempty(targetDoseIx)
                referenceDose = inf;
            else
                referenceDose = targetDoseInfo(targetDoseIx).refDose;
            end

            if ~isfinite(referenceDose)
                voiPrint = sprintf('%s%s', voiPrint, ...
                                   'Warning: target has no objective that penalizes underdosage, ');
            else
                StringReferenceDose = regexprep(num2str(round(referenceDose * 100) / 100), '\D', '_');
                VTarget95 = sum(doseInVoi >= 0.95 * referenceDose);
                VTreated95 = sum(doseCube(:) >= 0.95 * referenceDose);
                qi(runVoi).(['CI_' StringReferenceDose 'Gy']) = ...
                    VTarget95^2 / (numOfVoxels * VTreated95);

                qi(runVoi).(['HI_' StringReferenceDose 'Gy']) = ...
                    (DX(5) - DX(95)) / referenceDose * 100;
                qi(runVoi).referenceDose = referenceDose;
                qi(runVoi).COV_95 = VX(0.95 * referenceDose);
                qi(runVoi).COV_98 = VX(0.98 * referenceDose);
                qi(runVoi).COV_99 = VX(0.99 * referenceDose);
                qi(runVoi).COV1 = VX(referenceDose);

                voiPrint = sprintf(['%sCI = %6.4f, HI = %5.2f for ', ...
                                    'reference dose of %3.1f Gy\n%27s'], ...
                                   voiPrint, qi(runVoi).(['CI_' StringReferenceDose 'Gy']), ...
                                   qi(runVoi).(['HI_' StringReferenceDose 'Gy']), ...
                                   referenceDose, ' ');
                voiPrint = sprintf(['%sCOV95 = %6.2f%%, COV98 = %6.2f%%, ', ...
                                    'COV99 = %6.2f%%, COV1 = %6.2f%%\n'], ...
                                   voiPrint, 100 * qi(runVoi).COV_95, 100 * qi(runVoi).COV_98, ...
                                   100 * qi(runVoi).COV_99, 100 * qi(runVoi).COV1);
            end
        end

        matRad_cfg.dispInfo('%s\n', voiPrint);
    else
        matRad_cfg.dispInfo('%d %s - No dose information.', cst{runVoi, 1}, cst{runVoi, 2});
    end
end

listOfFields = fieldnames(qi);
for i = 1:size(cst, 1)
    indices = cst{i, 4}{1};
    doseInVoi = sort(doseCube(indices));
    if isempty(doseInVoi)
        for j = 1:numel(listOfFields)
            qi(i).(listOfFields{j}) = NaN;
        end
        qi(i).name = cst{i, 2};
    end
end

end
