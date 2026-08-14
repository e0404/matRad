function resultGUI = matRad_calcDeliveryMetrics(resultGUI, stf)
% matRad delivery metric calculation for aperture-based (VMAT) plans
%
% call
%   resultGUI = matRad_calcDeliveryMetrics(resultGUI, stf)
%
% input
%   resultGUI:  result struct carrying resultGUI.apertureInfo from
%               sequencing / direct aperture optimization
%   stf:        matRad steering struct (provides the FMO arc sectors used
%               for the leaf-travel direction statistics of continuous
%               aperture plans)
%
% output
%   resultGUI:  input struct with resultGUI.deliveryMetrics added::
%
%                 .planMU               total monitor units of the plan
%                 .planTime             total delivery time [s]
%                 .gantryAngle          gantry angle per DAO control point [deg]
%                 .time                 delivery time per DAO control point [s]
%                 .MURate               MU rate per DAO control point [MU/s]
%                 .gantryRotationSpeed  gantry speed per DAO control point [deg/s]
%                 .maxLeafSpeed         max leaf speed per DAO control point [mm/s]
%                 .leafSpeed            individual leaf speeds between control
%                                       points [mm/s] (signed for continuous
%                                       aperture delivery, absolute otherwise)
%                 .frac*                time fractions spent at the delivery
%                                       constraint limits
%
%               Continuous-aperture plans additionally get the leaf sweep
%               direction statistics .fracForward/.fracBackward per FMO
%               sector and their plan totals.
%
% References
%   [1] http://dx.doi.org/10.1118/1.4914863
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

if ~isfield(resultGUI, 'apertureInfo')
    matRad_cfg.dispError('Delivery metrics require resultGUI.apertureInfo from sequencing / direct aperture optimization!');
end
apertureInfo = resultGUI.apertureInfo;

if ~apertureInfo.runVMAT
    matRad_cfg.dispWarning('Delivery metrics are currently only computed for VMAT plans.');
    return
end

if ~isfield(apertureInfo, 'deliveryConstraints')
    matRad_cfg.dispError('apertureInfo carries no delivery constraints - run the VMAT sequencing first!');
end
constraints = apertureInfo.deliveryConstraints;

metrics = struct();
metrics.planMU   = 0;
metrics.planTime = 0;

% per-DAO-control-point quantities
numOfDAOBeams       = apertureInfo.totalNumOfShapes;
gantryAngle         = zeros(1, numOfDAOBeams);
time                = zeros(1, numOfDAOBeams);
MURate              = zeros(1, numOfDAOBeams);
gantryRotationSpeed = zeros(1, numOfDAOBeams);
maxLeafSpeed        = zeros(1, numOfDAOBeams);

l = 0;
for i = 1:numel(apertureInfo.beam)

    metrics.planMU   = metrics.planMU + apertureInfo.beam(i).shape(1).MU;
    metrics.planTime = metrics.planTime + apertureInfo.beam(i).time; % time until the next DAO control point

    if apertureInfo.beam(i).numOfShapes % only DAO control points carry their own delivery data
        l = l + 1;
        gantryAngle(l)         = apertureInfo.beam(i).gantryAngle;
        time(l)                = apertureInfo.beam(i).time;
        MURate(l)              = apertureInfo.beam(i).shape(1).MURate;
        gantryRotationSpeed(l) = apertureInfo.beam(i).gantryRot;
        maxLeafSpeed(l)        = apertureInfo.beam(i).maxLeafSpeed;
    end
end

metrics.gantryAngle         = gantryAngle;
metrics.time                = time;
metrics.MURate              = MURate;
metrics.gantryRotationSpeed = gantryRotationSpeed;
metrics.maxLeafSpeed        = maxLeafSpeed;

% individual leaf speeds between neighboring DAO control points, derived
% from the aperture vector
apertureInfoVec = apertureInfo.apertureVector;

% shorthand aliases
nShapes = apertureInfo.totalNumOfShapes;
nPairs  = apertureInfo.beam(1).numOfActiveLeafPairs;
nLP     = apertureInfo.totalNumOfLeafPairs;

leftLeafPos  = apertureInfoVec((1:nLP) + nShapes);
rightLeafPos = apertureInfoVec(1 + nLP + nShapes:nShapes + nLP * 2);
timeOptBorderAngles = apertureInfoVec((1 + nShapes + nLP * 2):end);

if apertureInfo.continuousAperture
    timeDoseBorderAngles = timeOptBorderAngles .* [apertureInfo.arc.beam([apertureInfo.arc.beam.isDAOBeam]).timeFactorCurrent]';

    leftLeafDiff = diff(reshape(leftLeafPos, nPairs, []), 1, 2);
    rightLeafDiff = diff(reshape(rightLeafPos, nPairs, []), 1, 2);

    isDAO = repmat([apertureInfo.arc.beam.isDAOBeam], nPairs, 1);
    leftLeafDiff  = reshape(leftLeafDiff(isDAO), nPairs, nShapes);
    rightLeafDiff = reshape(rightLeafDiff(isDAO), nPairs, nShapes);

    % signed leaf speeds: the sign encodes the travel direction
    metrics.leafSpeed = reshape([leftLeafDiff rightLeafDiff] ./ ...
                                repmat(timeDoseBorderAngles', nPairs, 2), 2 * nPairs * numel(timeDoseBorderAngles), 1);

    optAngles = [apertureInfo.beam([apertureInfo.arc.beam.isDAOBeam]).gantryAngle];
    metrics.leafSpeedGantryAngle = reshape(repmat(optAngles, nPairs, 2), 2 * nPairs * numel(timeDoseBorderAngles), 1);
else
    optInd = [apertureInfo.arc.beam.isDAOBeam];

    i = repelem(1:(nShapes - 1), 2);
    j = repelem(1:nShapes, 2);
    j(1) = [];
    j(end) = [];

    timeFactors = [apertureInfo.arc.beam(optInd).timeFactors]';
    timeFactors(1) = [];
    timeFactors(end) = [];

    timeFacMatrix = sparse(i, j, timeFactors, nShapes - 1, nShapes);
    timeBNOptAngles = timeFacMatrix * timeOptBorderAngles;

    metrics.leafSpeed = reshape([abs(diff(reshape(leftLeafPos, nPairs, nShapes), 1, 2)) ...
                                 abs(diff(reshape(rightLeafPos, nPairs, nShapes), 1, 2))] ./ ...
                                repmat(timeBNOptAngles', nPairs, 2), 2 * nPairs * numel(timeBNOptAngles), 1);
end

% time fractions spent at the delivery constraint limits (small relative
% tolerance against floating point comparison exactly at the limit)
tol = 1e-5;
totalTime = sum(time);
metrics.fracMaxMURate        = sum(time(MURate > constraints.monitorUnitRate(2) * (1 - tol))) ./ totalTime;
metrics.fracMinMURate        = sum(time(MURate < constraints.monitorUnitRate(1) * (1 + tol))) ./ totalTime;
metrics.fracMaxGantryRot     = sum(time(gantryRotationSpeed > constraints.gantryRotationSpeed(2) * (1 - tol))) ./ totalTime;
metrics.fracMaxLeafSpeed     = sum(time(maxLeafSpeed > constraints.leafSpeed(2) * (1 - tol))) ./ totalTime;
metrics.fracHalfMaxLeafSpeed = sum(time(maxLeafSpeed > constraints.leafSpeed(2) * (1 - tol) / 2)) ./ totalTime;

if apertureInfo.continuousAperture
    % leaf sweep direction statistics per FMO sector: the sweep direction
    % alternates between neighboring sectors
    if nargin < 2 || isempty(stf)
        matRad_cfg.dispError('The leaf sweep statistics of a continuous aperture plan require the stf!');
    end

    FMOBorders = [];
    for i = 1:numel(stf)
        if stf(i).arc.isFMOBeam
            FMOBorders = [FMOBorders, stf(i).arc.FMOAngleBorders(:)']; %#ok<AGROW>
        end
    end
    FMOBorders = unique(FMOBorders);

    forwardDir  = 1 - 2 * mod(1:(numel(FMOBorders) - 1), 2);
    numForward  = zeros(numel(forwardDir), 1);
    numBackward = zeros(numel(forwardDir), 1);
    for sector = 1:(numel(FMOBorders) - 1)
        inSector = FMOBorders(sector) <= metrics.leafSpeedGantryAngle & ...
                   metrics.leafSpeedGantryAngle <= FMOBorders(sector + 1);
        sectorSpeeds = metrics.leafSpeed(inSector);

        numForward(sector)  = nnz(sectorSpeeds * forwardDir(sector) >= 0);
        numBackward(sector) = nnz(sectorSpeeds * forwardDir(sector) < 0);
    end

    metrics.fracForward  = numForward ./ (numForward + numBackward);
    metrics.fracBackward = 1 - metrics.fracForward;
    metrics.totalFracForward  = mean(metrics.fracForward);
    metrics.totalFracBackward = 1 - metrics.totalFracForward;
end

resultGUI.deliveryMetrics = metrics;

end
