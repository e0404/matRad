function result = matRad_calcDeliveryMetrics(result, pln, stf)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad delivery metric calculation
%
% call
%   matRad_calcDeliveryMetrics(result,pln)
%
% input
%   result:             result struct from fluence optimization/sequencing
%   pln:                matRad plan meta information struct
%
% output
%   All plans: total MU
%   VMAT plans: total time, leaf speed, MU rate, and gantry rotation speed
%   distributions
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

apertureInfo = result.apertureInfo;

l = 0;
if apertureInfo.runVMAT

    machine = matRad_loadMachine(pln);

    apertureInfo.planMU     = 0;
    apertureInfo.planTime   = 0;

    % All of these are vectors
    % Each entry corresponds to a beam angle
    % Later, we will convert these to histograms, find max, mean, min, etc.
    gantryRot = zeros(1, result.apertureInfo.totalNumOfShapes);
    MURate = gantryRot;
    times = gantryRot;
    angles = gantryRot;
    maxLeafSpeed = gantryRot;

    for i = 1:size(apertureInfo.beam, 2)

        apertureInfo.planMU     = apertureInfo.planMU + apertureInfo.beam(i).shape(1).MU;
        apertureInfo.planTime   = apertureInfo.planTime + apertureInfo.beam(i).time; % time until next optimized beam

        if apertureInfo.beam(i).numOfShapes % only optimized beams have their time in the data struct
            l = l + 1;
            gantryRot(l) = apertureInfo.beam(i).gantryRot;
            MURate(l) = apertureInfo.beam(i).shape(1).MURate * 60;
            times(l) = apertureInfo.beam(i).time;
            maxLeafSpeed(l) = apertureInfo.beam(i).maxLeafSpeed / 10;
            angles(l) = apertureInfo.beam(i).gantryAngle;
        end
    end

    apertureInfoVec = apertureInfo.apertureVector;

    % shorthand aliases
    nShapes = apertureInfo.totalNumOfShapes;
    nPairs  = apertureInfo.beam(1).numOfActiveLeafPairs;
    nLP     = apertureInfo.totalNumOfLeafPairs;

    leftLeafPos  = apertureInfoVec((1:nLP) + nShapes);
    rightLeafPos = apertureInfoVec(1 + nLP + nShapes:nShapes + nLP * 2);
    timeOptBorderAngles = apertureInfoVec((1 + nShapes + nLP * 2):end);

    if apertureInfo.continuousAperture
        timeDoseBorderAngles = timeOptBorderAngles .* [apertureInfo.arc.beam([apertureInfo.arc.beam.DAOBeam]).timeFacCurr]';

        leftLeafDiff = diff(reshape(leftLeafPos, nPairs, []), 1, 2);
        rightLeafDiff = diff(reshape(rightLeafPos, nPairs, []), 1, 2);

        isDAO = repmat([apertureInfo.arc.beam.DAOBeam], nPairs, 1);
        leftLeafDiff  = reshape(leftLeafDiff(isDAO), nPairs, nShapes);
        rightLeafDiff = reshape(rightLeafDiff(isDAO), nPairs, nShapes);

        lfspd = reshape([leftLeafDiff rightLeafDiff] ./ ...
                        repmat(timeDoseBorderAngles', nPairs, 2), 2 * nPairs * numel(timeDoseBorderAngles), 1);

        optAngles = [apertureInfo.beam([apertureInfo.arc.beam.DAOBeam]).gantryAngle];
        optAnglesMat = reshape(repmat(optAngles, nPairs, 2), 2 * nPairs * numel(timeDoseBorderAngles), 1);
    else
        optInd = [apertureInfo.arc.beam.DAOBeam];

        i = repelem(1:(nShapes - 1), 2);
        j = repelem(1:nShapes, 2);
        j(1) = [];
        j(end) = [];

        timeFac = [apertureInfo.arc.beam(optInd).timeFac]';
        timeFac(1) = [];
        timeFac(end) = [];

        timeFacMatrix = sparse(i, j, timeFac, nShapes - 1, nShapes);
        timeBNOptAngles = timeFacMatrix * timeOptBorderAngles;

        lfspd = reshape([abs(diff(reshape(leftLeafPos, nPairs, nShapes), 1, 2)) ...
                         abs(diff(reshape(rightLeafPos, nPairs, nShapes), 1, 2))] ./ ...
                        repmat(timeBNOptAngles', nPairs, 2), 2 * nPairs * numel(timeBNOptAngles), 1);
    end

    if apertureInfo.continuousAperture

        % FMOBorders = zeros(1,2*numel(pln.propStf.FMOGantryAngles));
        counter = 1;
        for i = 1:numel(stf)
            if stf(i).arc.FMOBeam
                FMOBorders(counter) = stf(i).arc.FMOAngleBorders(1);
                FMOBorders(counter + 1) = stf(i).arc.FMOAngleBorders(2);
                counter = counter + 2;
            else
                continue
            end
        end
        FMOBorders = unique(FMOBorders);
        forwardDir = 1 - 2 * mod(1:(numel(FMOBorders) - 1), 2);
        numForward = zeros(numel(forwardDir), 1);
        numBackward = zeros(numel(forwardDir), 1);
        timeInInit = zeros(numel(forwardDir), 1);

        plot(optAnglesMat, lfspd, '.');
        hold on;
        counter = 1;
        for border = FMOBorders
            plot([border border], [-machine.constraints.leafSpeed(2) machine.constraints.leafSpeed(2)], 'r-');

            if border < FMOBorders(end)
                curr_lfspd = lfspd(FMOBorders(counter) <= optAnglesMat & optAnglesMat <= FMOBorders(counter + 1));

                numForward(counter) = nnz(curr_lfspd * forwardDir(counter) >= 0);
                numBackward(counter) = nnz(curr_lfspd * forwardDir(counter) < 0);
                timeInInit(counter) = sum(times(FMOBorders(counter) <= angles & angles <= FMOBorders(counter + 1)));

                counter = counter + 1;
            end
        end

        figure;
        plot([min(FMOBorders) - 5 max(FMOBorders) + 5], [0 0], 'k--');
        xlim([min(FMOBorders) - 5 max(FMOBorders) + 5]);
        ylim([-machine.constraints.leafSpeed(2) - 5 machine.constraints.leafSpeed(2) + 5]);
        xlabel('gantry angle (^\circ)');
        ylabel('leaf speed (cm/s)');

        figure;
        plot(optAngles, gantryRot, '.');
        xlim([min(FMOBorders) - 5 max(FMOBorders) + 5]);
        ylim([0 machine.constraints.gantryRotationSpeed(2) + 1]);
        xlabel('gantry angle (^\circ)');
        ylabel('gantry rotation speed (^\circ/s)');

        figure;
        plot(optAngles, MURate, '.');
        xlim([min(FMOBorders) - 5 max(FMOBorders) + 5]);
        ylim([0 60 * machine.constraints.monitorUnitRate(2) + 5]);
        xlabel('gantry angle (^\circ)');
        ylabel('MU rate (MU/min)');

        apertureInfo.fracMaxMURate = sum(times(MURate > 60 * machine.constraints.monitorUnitRate(2) * (1 - 1e-5))) ./ sum(times);
        apertureInfo.fracMinMURate = sum(times(MURate < 60 * machine.constraints.monitorUnitRate(1) * (1 + 1e-5))) ./ sum(times);
        apertureInfo.fracMaxGantryRot = sum(times(gantryRot > machine.constraints.gantryRotationSpeed(2) * (1 - 1e-5))) ./ sum(times);
        apertureInfo.fracMaxLeafSpeed = sum(times(maxLeafSpeed > machine.constraints.leafSpeed(2) / 10 * (1 - 1e-5))) ./ sum(times);
        apertureInfo.fracHalfMaxLeafSpeed = sum(times(maxLeafSpeed > machine.constraints.leafSpeed(2) / 10 * (1 - 1e-5) / 2)) ./ sum(times);

        apertureInfo.fracForward = numForward ./ (numForward + numBackward);
        apertureInfo.fracBackward = 1 - apertureInfo.fracForward;
        apertureInfo.totalFracForward = mean(apertureInfo.fracForward);
        % apertureInfo.totalFracForward = sum(apertureInfo.fracForward.*timeInInit)./sum(timeInInit);
        apertureInfo.totalFracBackward = 1 - apertureInfo.totalFracForward;
    end
    %}

end

result.apertureInfo = apertureInfo;
