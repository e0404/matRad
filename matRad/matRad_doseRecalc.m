function recalc = matRad_doseRecalc(cst, pln, recalc, ct, apertureInfo, calcDoseDirect, dij)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to recalculate dose on a finer or equal angular
% resolution, either by interpolating aperture shapes or by reusing the
% nearest existing Dij column.
%
% call
%   recalc = matRad_doseRecalc(cst,pln,recalc,ct,apertureInfo)
%   recalc = matRad_doseRecalc(cst,pln,recalc,ct,apertureInfo,calcDoseDirect)
%   recalc = matRad_doseRecalc(cst,pln,recalc,ct,apertureInfo,calcDoseDirect,dij)
%
% input
%   cst, ct:          patient data
%   pln:              original optimisation plan (anchor angles + propOpt)
%   recalc:           recalc options struct (interpNew, dijNew,
%                     continuousAperture, pln with recalc spacing, ...)
%   apertureInfo:     aperture info from the optimisation result
%   calcDoseDirect:   (optional, default true) use direct dose calc
%   dij:              (optional) Dij for back-projection when ~calcDoseDirect
%
% output
%   recalc:  updated struct with stf, apertureInfo, and resultGUI
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
    calcDoseDirect = true;
end

recalc.apertureInfo = apertureInfo;

% Old fine angles are the DAO beams from the original optimisation.
% They live in apertureInfo, so we do not need pln.propStf for this.
oldFineAngles = [apertureInfo.beam.gantryAngle];

% If not interpolating new apertures, force every fine angle to be a DAO
% control point by equalising the two spacing parameters.
if ~recalc.interpNew
    recalc.pln.propStf.maxDAOGantryAngleSpacing = recalc.pln.propStf.maxGantryAngleSpacing;
end

% Generate (or load cached) stf at the recalc angular resolution.
% recalc.pln.propStf.gantryAngles contains arc anchor points; the fine
% angle grid is computed internally by matRad_StfGeneratorPhotonVMAT.
fname = sprintf('%.1f deg.mat', recalc.pln.propStf.maxGantryAngleSpacing);
if exist(fname, 'file')
    load(fname, 'stf');
else
    stf = matRad_generateStf(ct, cst, recalc.pln);
    save(fname, 'stf');
end
recalc.stf = stf;
clear stf;

% -----------------------------------------------------------------------
% Handle angles that fall exactly equidistant between two old DAO angles.
% These need duplicate stf entries so that both neighbouring Dij columns
% can be used (one per side).
% -----------------------------------------------------------------------
if ~recalc.interpNew || ~recalc.dijNew

    newFineAngles = [recalc.stf.gantryAngle];

    duplicate = false(size(newFineAngles));
    for i = 1:numel(newFineAngles)
        diffs = abs(newFineAngles(i) - oldFineAngles);
        duplicate(i) = sum(diffs == min(diffs)) > 1;
    end

    tempStf = recalc.stf;
    [tempStf(:).copyInd] = deal([]);
    [tempStf(:).stfCorr]  = deal([]);

    j = 1;
    for i = 1:numel(newFineAngles)
        if duplicate(i)
            % Left copy: pretend this beam sits at the previous beam angle
            tempStf(j)            = recalc.stf(i);
            tempStf(j).gantryAngle = recalc.stf(i - 1).gantryAngle;
            tempStf(j).copyInd    = 1;
            tempStf(j).stfCorr    = false;
            j = j + 1;

            % Right copy: pretend this beam sits at the next beam angle
            tempStf(j)            = recalc.stf(i);
            tempStf(j).gantryAngle = recalc.stf(i + 1).gantryAngle;
            tempStf(j).copyInd    = 2;
            tempStf(j).stfCorr    = false;
        else
            tempStf(j)         = recalc.stf(i);
            tempStf(j).copyInd = [];
            tempStf(j).stfCorr = true;
        end
        j = j + 1;
    end
    recalc.stf = tempStf;

    % -------------------------------------------------------------------
    % Dij reuse: redirect each new beam to the nearest old beam so that
    % an existing Dij column can be recycled without recomputation.
    % -------------------------------------------------------------------
    tempStf = recalc.stf;
    for i = 1:numel(recalc.stf)
        diffs      = abs(recalc.stf(i).gantryAngle - oldFineAngles);
        minDiff    = min(diffs);
        nearIdx    = find(diffs == minDiff);   % 1 or 2 indices into oldFineAngles

        % Find where those old angles appear in the (post-duplicate) stf
        newAngles  = [tempStf.gantryAngle];
        minInd1    = find(newAngles == oldFineAngles(nearIdx(1)),   1);
        minInd2    = find(newAngles == oldFineAngles(nearIdx(end)), 1);

        if ~recalc.dijNew
            % Replace this beam with the nearest existing beam entirely
            if isempty(recalc.stf(i).copyInd) || recalc.stf(i).copyInd == 1
                recalc.stf(i) = tempStf(minInd1);
            elseif recalc.stf(i).copyInd == 2
                recalc.stf(i) = tempStf(minInd2);
            end
        elseif ~recalc.interpNew
            % Keep the beam but correct its angle if it is equidistant
            if numel(nearIdx) > 1
                recalc.stf(i).gantryAngle = tempStf(i).gantryAngle;
            end
        end
    end
end

recalc = matRad_recalcApertureInfo(recalc, recalc.apertureInfo);

recalc.apertureInfo.propVMAT.continuousAperture = recalc.continuousAperture;
recalc.apertureInfo = matRad_daoVec2ApertureInfo_VMATrecalcDynamic( ...
                                                                   recalc.apertureInfo, recalc.apertureInfo.apertureVector);

if calcDoseDirect
    clear global;
    recalc.resultGUI = matRad_calcDoseDirect(ct, recalc.stf, recalc.pln, cst, ...
                                             recalc.apertureInfo.bixelWeights);
else
    recalc.resultGUI.w = recalc.apertureInfo.bixelWeights;

    options.numOfScenarios = 1;
    options.bioOpt         = 'none';
    dij.scaleFactor        = apertureInfo.weightToMU ./ dij.weightToMU;
    d = matRad_backProjection(recalc.resultGUI.w, dij, options);
    recalc.resultGUI.physicalDose = reshape(d{1}, dij.dimensions);
end
