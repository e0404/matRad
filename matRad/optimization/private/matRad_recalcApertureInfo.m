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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stf = recalc.stf;

% accept aperture info saved by an older matRad version
apertureInfoOld = matRad_upgradeApertureInfo(apertureInfoOld);

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

    % Build a leaf trajectory at dose-sector borders. Continuous plans
    % carry the actual initial/final leaf positions there. For a discrete
    % source plan, extend the first/last aperture to the arc boundaries.
    if matRad_getFieldOrDefault(apertureInfoOld, 'continuousAperture', false)
        oldLeftFinal = arrayfun(@(b) b.shape(1).leftLeafPosFinal, apertureInfoOld.beam, ...
                                'UniformOutput', false);
        oldRightFinal = arrayfun(@(b) b.shape(1).rightLeafPosFinal, apertureInfoOld.beam, ...
                                 'UniformOutput', false);
        oldBorderAngles = [apertureInfoOld.arc.beam(1).doseAngleBorders(1), ...
                           arrayfun(@(b) b.doseAngleBorders(2), apertureInfoOld.arc.beam)];
        oldLeftBorderPoss = [apertureInfoOld.beam(1).shape(1).leftLeafPosInitial, ...
                             oldLeftFinal{:}];
        oldRightBorderPoss = [apertureInfoOld.beam(1).shape(1).rightLeafPosInitial, ...
                              oldRightFinal{:}];
    else
        oldBorderAngles = [apertureInfoOld.arc.beam(1).doseAngleBorders(1), oldGantryAngles, ...
                           apertureInfoOld.arc.beam(end).doseAngleBorders(2)];
        oldLeftBorderPoss = [oldLeftLeafPoss(:, 1), oldLeftLeafPoss, oldLeftLeafPoss(:, end)];
        oldRightBorderPoss = [oldRightLeafPoss(:, 1), oldRightLeafPoss, oldRightLeafPoss(:, end)];
    end

    % The new beam centres are not bracketed by the old beam centres: the
    % outermost new centres sit between the arc boundary and the first/last
    % old centre. Anchor the centre trajectory at the arc boundaries too, so
    % those beams interpolate instead of extrapolating to NaN.
    oldCentreAngles = [oldBorderAngles(1), oldGantryAngles, oldBorderAngles(end)];
    oldLeftCentrePoss = [oldLeftBorderPoss(:, 1), oldLeftLeafPoss, oldLeftBorderPoss(:, end)];
    oldRightCentrePoss = [oldRightBorderPoss(:, 1), oldRightLeafPoss, oldRightBorderPoss(:, end)];
end

% MLC parameters:
numOfMLCLeafPairs = 80;

% initializing variables
totalNumOfShapes = numel(stf);
for i = 1:numel(apertureInfoOld.beam)
    oldSectorBounds = sort(apertureInfoOld.arc.beam(i).doseAngleBorders);
    newInd = (oldSectorBounds(1) <= [stf.gantryAngle] & ...
              [stf.gantryAngle] <= oldSectorBounds(2)) .* ...
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
        apertureInfoNew.beam(j).limLeft = geometry.limLeft;
        apertureInfoNew.beam(j).limRight = geometry.limRight;
        apertureInfoNew.beam(j).bixelIndMap = geometry.bixelIndMap;
        apertureInfoNew.beam(j).posOfCornerBixel = geometry.posOfCornerBixel;
        apertureInfoNew.beam(j).MLCWindow = geometry.MLCWindow;
        apertureInfoNew.beam(j).bixOffset = 1 + (j - 1) * dimZ;
        apertureInfoNew.beam(j).shape(1).vectorOffset = totalNumOfShapes + 1 + (j - 1) * dimZ;

        % inherit from old beam
        apertureInfoNew.arc.beam(j).leafDir = apertureInfoOld.arc.beam(i).leafDir;

        % specific to new beam
        apertureInfoNew.beam(j).gantryAngle = stf(j).gantryAngle;
        apertureInfoNew.arc.beam(j).doseAngleBorders = stf(j).arc.doseAngleBorders;
        apertureInfoNew.arc.beam(j).doseAngleBorderCentreDiff = stf(j).arc.doseAngleBorderCentreDiff;
        apertureInfoNew.arc.beam(j).doseAngleBordersDiff = stf(j).arc.doseAngleBordersDiff;
        apertureInfoNew.arc.beam(j).lastDAOBeamIx = stf(j).arc.lastDAOBeamIx;
        apertureInfoNew.arc.beam(j).nextDAOBeamIx = stf(j).arc.nextDAOBeamIx;

        newSectorBounds = sort(apertureInfoNew.arc.beam(j).doseAngleBorders);
        overlap = matRad_intervalOverlap(newSectorBounds, oldSectorBounds);
        amountOfOldSpeed = overlap ./ apertureInfoNew.arc.beam(j).doseAngleBordersDiff;
        amountOfOldWeight = overlap ./ apertureInfoOld.arc.beam(i).doseAngleBordersDiff;

        totalAmountOfOldWeight = totalAmountOfOldWeight + amountOfOldWeight;

        initialHalfBounds = sort([apertureInfoNew.arc.beam(j).doseAngleBorders(1), ...
                                  apertureInfoNew.beam(j).gantryAngle]);
        finalHalfBounds = sort([apertureInfoNew.beam(j).gantryAngle, ...
                                apertureInfoNew.arc.beam(j).doseAngleBorders(2)]);
        amountOfOldWeight_I = matRad_intervalOverlap(initialHalfBounds, oldSectorBounds) ./ ...
                              apertureInfoOld.arc.beam(i).doseAngleBordersDiff;
        amountOfOldWeight_F = matRad_intervalOverlap(finalHalfBounds, oldSectorBounds) ./ ...
                              apertureInfoOld.arc.beam(i).doseAngleBordersDiff;

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
            apertureInfoNew.arc.beam(j).doseAngleBordersDiff;

        % apertureInfoNew.beam(j).shape(1).jacobiScale = apertureInfoOld.beam(i).shape(1).jacobiScale;
        apertureInfoNew.jacobiScale(j) = 1;
        apertureInfoNew.beam(j).shape(1).jacobiScale = 1;

        if recalc.interpNew
            % interpolate new apertures now so that weights are not
            % overwritten
            beamAngle = apertureInfoNew.beam(j).gantryAngle;
            doseAngleBorders = apertureInfoNew.arc.beam(j).doseAngleBorders;

            leftLeafPos = matRad_interpLeafTrajectory(oldCentreAngles, oldLeftCentrePoss, beamAngle);
            rightLeafPos = matRad_interpLeafTrajectory(oldCentreAngles, oldRightCentrePoss, beamAngle);
            leftLeafPosInitial = matRad_interpLeafTrajectory(oldBorderAngles, oldLeftBorderPoss, doseAngleBorders(1));
            rightLeafPosInitial = matRad_interpLeafTrajectory(oldBorderAngles, oldRightBorderPoss, doseAngleBorders(1));
            leftLeafPosFinal = matRad_interpLeafTrajectory(oldBorderAngles, oldLeftBorderPoss, doseAngleBorders(2));
            rightLeafPosFinal = matRad_interpLeafTrajectory(oldBorderAngles, oldRightBorderPoss, doseAngleBorders(2));

            apertureInfoNew.beam(j).shape(1).leftLeafPos = leftLeafPos;
            apertureInfoNew.beam(j).shape(1).rightLeafPos = rightLeafPos;
            apertureInfoNew.beam(j).shape(1).leftLeafPosInitial = leftLeafPosInitial;
            apertureInfoNew.beam(j).shape(1).rightLeafPosInitial = rightLeafPosInitial;
            apertureInfoNew.beam(j).shape(1).leftLeafPosFinal = leftLeafPosFinal;
            apertureInfoNew.beam(j).shape(1).rightLeafPosFinal = rightLeafPosFinal;
        else
            apertureInfoNew.beam(j).shape(1).leftLeafPos = apertureInfoOld.beam(i).shape(1).leftLeafPos;
            apertureInfoNew.beam(j).shape(1).rightLeafPos = apertureInfoOld.beam(i).shape(1).rightLeafPos;
        end

        % all beams are now "optimized" to prevent their weights from being
        % overwritten
        % optAngleBorders becomes doseAngleBorders
        apertureInfoNew.beam(j).numOfShapes = 1;
        apertureInfoNew.arc.beam(j).isDAOBeam = true;
        apertureInfoNew.arc.beam(j).DAOAngleBorders = stf(j).arc.doseAngleBorders;
        apertureInfoNew.arc.beam(j).DAOAngleBorderCentreDiff = stf(j).arc.doseAngleBorderCentreDiff;
        apertureInfoNew.arc.beam(j).DAOAngleBordersDiff = stf(j).arc.doseAngleBordersDiff;
        apertureInfoNew.arc.beam(j).timeFactorCurrent = ...
            apertureInfoNew.arc.beam(j).doseAngleBordersDiff ./ apertureInfoNew.arc.beam(j).DAOAngleBordersDiff; % = 1

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

end

function leafPos = matRad_interpLeafTrajectory(angles, positions, queryAngle)
% Interpolate a leaf trajectory independent of the arc direction. Duplicate
% angle samples can occur for step-and-shoot endpoints and are collapsed.

[angles, sortIx] = sort(angles);
positions = positions(:, sortIx);
[angles, uniqueIx] = unique(angles, 'stable');
positions = positions(:, uniqueIx);

leafPos = interp1(angles', positions', queryAngle, 'linear')';

end

function overlap = matRad_intervalOverlap(firstBounds, secondBounds)
% Length of the overlap between two direction-independent angle intervals.

overlap = max(0, min(firstBounds(2), secondBounds(2)) - ...
              max(firstBounds(1), secondBounds(1)));

end
