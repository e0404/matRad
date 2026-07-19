function recalc = matRad_recalcApertureInfo(recalc, apertureInfoOld)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to apertures for a different dose resolution
%
% call
%   recalc = matRad_recalcApertureInfo(recalc,apertureInfo)
%
% input
%   recalc:
%   apertureInfo:
%
% output
%   recalc:
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

stf = recalc.stf;

apertureInfoNew = apertureInfoOld;
apertureInfoNew = rmfield(apertureInfoNew, 'beam');

apertureInfoNew.totalNumOfBixels = sum([stf(:).totalNumOfBixels]);

shapeInd = 1;

if recalc.interpNew
    oldGantryAngles = zeros(1, numel(apertureInfoOld.beam));
    oldLeftLeafPoss = zeros(apertureInfoOld.beam(1).numOfActiveLeafPairs, numel(apertureInfoOld.beam));
    oldRightLeafPoss = zeros(apertureInfoOld.beam(1).numOfActiveLeafPairs, numel(apertureInfoOld.beam));
    for i = 1:numel(apertureInfoOld.beam)
        oldGantryAngles(i) = apertureInfoOld.beam(i).gantryAngle;
        oldLeftLeafPoss(:, i) = apertureInfoOld.beam(i).shape(1).leftLeafPos;
        oldRightLeafPoss(:, i) = apertureInfoOld.beam(i).shape(1).rightLeafPos;
    end
end

% MLC parameters:
numOfMLCLeafPairs = 80;

% initializing variables
totalNumOfShapes = numel(stf);
for i = 1:numel(apertureInfoOld.beam)
    newInd = (apertureInfoOld.propVMAT.beam(i).doseAngleBorders(1) <= [stf.gantryAngle] & ...
              [stf.gantryAngle] <= apertureInfoOld.propVMAT.beam(i).doseAngleBorders(2)) .* ...
              (1:numel([stf.gantryAngle]));
    newInd(newInd == 0) = [];

    totalAmountOfOldWeight = 0;

    for j = newInd
        % derive the MLC geometry for this beam; the bixel numbering
        % restarts per beam index (equal ray count per beam)
        geometry = matRad_PhotonSequencerAbstract.getMLCGeometry(stf(j), numOfMLCLeafPairs, (j - 1) * stf(1).numOfRays);
        dimZ = geometry.numOfActiveLeafPairs;

        % save data for each beam
        apertureInfoNew.beam(j).numOfActiveLeafPairs = dimZ;
        apertureInfoNew.beam(j).leafPairPos = geometry.leafPairPos;
        apertureInfoNew.beam(j).isActiveLeafPair = geometry.isActiveLeafPair;
        apertureInfoNew.beam(j).centralLeafPair = geometry.centralLeafPair;
        apertureInfoNew.beam(j).lim_l = geometry.lim_l;
        apertureInfoNew.beam(j).lim_r = geometry.lim_r;
        apertureInfoNew.beam(j).bixelIndMap = geometry.bixelIndMap;
        apertureInfoNew.beam(j).posOfCornerBixel = geometry.posOfCornerBixel;
        apertureInfoNew.beam(j).MLCWindow = geometry.MLCWindow;
        apertureInfoNew.beam(j).bixOffset = 1 + (j - 1) * dimZ;
        apertureInfoNew.beam(j).shape(1).vectorOffset = totalNumOfShapes + 1 + (j - 1) * dimZ;

        % inherit from old beam
        apertureInfoNew.propVMAT.beam(j).leafDir = apertureInfoOld.propVMAT.beam(i).leafDir;

        % specific to new beam
        apertureInfoNew.beam(j).gantryAngle = stf(j).gantryAngle;
        apertureInfoNew.propVMAT.beam(j).doseAngleBorders = stf(j).propVMAT.doseAngleBorders;
        apertureInfoNew.propVMAT.beam(j).doseAngleBorderCentreDiff = stf(j).propVMAT.doseAngleBorderCentreDiff;
        apertureInfoNew.propVMAT.beam(j).doseAngleBordersDiff = stf(j).propVMAT.doseAngleBordersDiff;
        apertureInfoNew.propVMAT.beam(j).lastDAOIndex = stf(j).propVMAT.lastDAOIndex;
        apertureInfoNew.propVMAT.beam(j).nextDAOIndex = stf(j).propVMAT.lastDAOIndex;

        amountOfOldSpeed = (min(apertureInfoNew.propVMAT.beam(j).doseAngleBorders(2), ...
                                apertureInfoOld.propVMAT.beam(i).doseAngleBorders(2)) - ...
                            max(apertureInfoNew.propVMAT.beam(j).doseAngleBorders(1), ...
                                apertureInfoOld.propVMAT.beam(i).doseAngleBorders(1))) ./ ...
            apertureInfoNew.propVMAT.beam(j).doseAngleBordersDiff;
        amountOfOldWeight = (min(apertureInfoNew.propVMAT.beam(j).doseAngleBorders(2), ...
                                 apertureInfoOld.propVMAT.beam(i).doseAngleBorders(2)) - ...
                             max(apertureInfoNew.propVMAT.beam(j).doseAngleBorders(1), ...
                                 apertureInfoOld.propVMAT.beam(i).doseAngleBorders(1))) ./ ...
            apertureInfoOld.propVMAT.beam(i).doseAngleBordersDiff;

        totalAmountOfOldWeight = totalAmountOfOldWeight + amountOfOldWeight;

        amountOfOldWeight_I = (min(apertureInfoNew.beam(j).gantryAngle, ...
                                   apertureInfoOld.propVMAT.beam(i).doseAngleBorders(2)) - ...
                               max(apertureInfoNew.propVMAT.beam(j).doseAngleBorders(1), ...
                                   apertureInfoOld.propVMAT.beam(i).doseAngleBorders(1))) ./ ...
                                     apertureInfoOld.propVMAT.beam(i).doseAngleBordersDiff;
        amountOfOldWeight_F = (min(apertureInfoNew.propVMAT.beam(j).doseAngleBorders(2), ...
                                   apertureInfoOld.propVMAT.beam(i).doseAngleBorders(2)) - ...
                               max(apertureInfoNew.beam(j).gantryAngle, ...
                                   apertureInfoOld.propVMAT.beam(i).doseAngleBorders(1))) ./ ...
                                     apertureInfoOld.propVMAT.beam(i).doseAngleBordersDiff;

        if ~isfield(apertureInfoNew.beam(j), 'gantryRot') || isempty(apertureInfoNew.beam(j).gantryRot)
            apertureInfoNew.beam(j).gantryRot = 0;
            apertureInfoNew.beam(j).shape(1).weight = 0;
            apertureInfoNew.beam(j).shape(1).weight_I = 0;
            apertureInfoNew.beam(j).shape(1).weight_F = 0;
        end
        apertureInfoNew.beam(j).gantryRot = amountOfOldSpeed * apertureInfoOld.beam(i).gantryRot + apertureInfoNew.beam(j).gantryRot;

        % recalculate weight, MU
        apertureInfoNew.beam(j).shape(1).weight = apertureInfoNew.beam(j).shape(1).weight + ...
            amountOfOldWeight * apertureInfoOld.beam(i).shape(1).weight;
        apertureInfoNew.beam(j).shape(1).weight_I = apertureInfoNew.beam(j).shape(1).weight_I + ...
            amountOfOldWeight_I * apertureInfoOld.beam(i).shape(1).weight;
        apertureInfoNew.beam(j).shape(1).weight_F = apertureInfoNew.beam(j).shape(1).weight_F + ...
            amountOfOldWeight_F * apertureInfoOld.beam(i).shape(1).weight;
        apertureInfoNew.beam(j).MU = apertureInfoNew.beam(j).shape(1).weight .* ...
            apertureInfoNew.weightToMU;

        apertureInfoNew.beam(j).MURate = apertureInfoNew.beam(j).MU .* apertureInfoNew.beam(j).gantryRot ./ ...
            apertureInfoNew.propVMAT.beam(j).doseAngleBordersDiff;

        % apertureInfoNew.beam(j).shape(1).jacobiScale = apertureInfoOld.beam(i).shape(1).jacobiScale;
        apertureInfoNew.jacobiScale(j) = 1;
        apertureInfoNew.beam(j).shape(1).jacobiScale = 1;

        if recalc.interpNew
            % interpolate new apertures now so that weights are not
            % overwritten
            apertureInfoNew.beam(j).shape(1).leftLeafPos = ...
                (interp1(oldGantryAngles', oldLeftLeafPoss', apertureInfoNew.beam(j).gantryAngle))';
            apertureInfoNew.beam(j).shape(1).rightLeafPos = ...
                (interp1(oldGantryAngles', oldRightLeafPoss', apertureInfoNew.beam(j).gantryAngle))';

            apertureInfoNew.beam(j).shape(1).leftLeafPos_I = ...
                (interp1(oldGantryAngles', oldLeftLeafPoss', apertureInfoNew.propVMAT.beam(j).doseAngleBorders(1)))';
            apertureInfoNew.beam(j).shape(1).rightLeafPos_I = ...
                (interp1(oldGantryAngles', oldRightLeafPoss', apertureInfoNew.propVMAT.beam(j).doseAngleBorders(1)))';

            apertureInfoNew.beam(j).shape(1).leftLeafPos_F = ...
                (interp1(oldGantryAngles', oldLeftLeafPoss', apertureInfoNew.propVMAT.beam(j).doseAngleBorders(2)))';
            apertureInfoNew.beam(j).shape(1).rightLeafPos_F = ...
                (interp1(oldGantryAngles', oldRightLeafPoss', apertureInfoNew.propVMAT.beam(j).doseAngleBorders(2)))';
        else
            apertureInfoNew.beam(j).shape(1).leftLeafPos = apertureInfoOld.beam(i).shape(1).leftLeafPos;
            apertureInfoNew.beam(j).shape(1).rightLeafPos = apertureInfoOld.beam(i).shape(1).rightLeafPos;
        end

        % all beams are now "optimized" to prevent their weights from being
        % overwritten
        % optAngleBorders becomes doseAngleBorders
        apertureInfoNew.beam(j).numOfShapes = 1;
        apertureInfoNew.propVMAT.beam(j).DAOBeam = true;
        apertureInfoNew.propVMAT.beam(j).DAOAngleBorders = stf(j).propVMAT.doseAngleBorders;
        apertureInfoNew.propVMAT.beam(j).DAOAngleBorderCentreDiff = stf(j).propVMAT.doseAngleBorderCentreDiff;
        apertureInfoNew.propVMAT.beam(j).DAOAngleBordersDiff = stf(j).propVMAT.doseAngleBordersDiff;
        apertureInfoNew.propVMAT.beam(j).timeFacCurr = ...
            apertureInfoNew.propVMAT.beam(j).doseAngleBordersDiff ./ apertureInfoNew.propVMAT.beam(j).DAOAngleBordersDiff; % = 1

        apertureInfoNew.apertureVector(shapeInd) = apertureInfoNew.beam(j).shape(1).weight;
        shapeInd = shapeInd + 1;
    end

end

apertureInfoNew.totalNumOfShapes        = sum([apertureInfoNew.beam.numOfShapes]);
apertureInfoNew.totalNumOfLeafPairs     = sum([apertureInfoNew.beam.numOfShapes] * [apertureInfoNew.beam.numOfActiveLeafPairs]');
apertureInfoNew.doseTotalNumOfLeafPairs = sum([apertureInfoNew.beam(:).numOfActiveLeafPairs]);
apertureInfoNew.totalNumOfOptBixels     = apertureInfoNew.totalNumOfBixels;

% The rebuilt beams always use the discrete (single leaf position set)
% vector layout - the sweep mode for the recalculated delivery is applied
% separately in matRad_recalcApertureBixelWeights.
apertureInfoNew.continuousAperture = false;

% recalc apertureVector
[apertureInfoNew.apertureVector, apertureInfoNew.mappingMx, apertureInfoNew.limMx] = ...
    matRad_OptimizationProblemVMAT.matRad_daoApertureInfo2Vec(apertureInfoNew);

recalc.apertureInfo = apertureInfoNew;
recalc.stf = stf;
