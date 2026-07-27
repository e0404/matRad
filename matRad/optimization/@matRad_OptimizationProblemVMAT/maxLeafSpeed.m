function apertureInfo = maxLeafSpeed(apertureInfo)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad calculation of maximum leaf speed
%
% call
%   apertureInfo = matRad_OptimizationProblemVMAT.maxLeafSpeed(apertureInfo)
%
% input
%   apertureInfo:   aperture info struct
%
% output
%   apertureInfo:   aperture info struct
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

apertureInfoVec = apertureInfo.apertureVector;

% values of time differences of optimized gantry angles
timeDAOBorderAngles = apertureInfoVec(1 + (apertureInfo.totalNumOfShapes + apertureInfo.totalNumOfLeafPairs * 2):end);

% find values of leaf speeds of optimized gantry angles
if apertureInfo.continuousAperture
    % Using the dynamic fluence calculation, we have the leaf positions in
    % the vector be the leaf positions at the borders of the Dij arcs (for optimized angles only).
    % Therefore we must also use the times between the borders of the Dij
    % arc (for optimized angles only).
    timeFactors = [apertureInfo.arc.beam.timeFactors]';
    deleteInd = timeFactors == 0;
    timeFactors(deleteInd) = [];

    i = [apertureInfo.arc.beam.timeFactorIx]';
    i(deleteInd) = [];

    j = repelem(1:apertureInfo.totalNumOfShapes, 1, 3);
    j(deleteInd) = [];

    timeFacMatrix = sparse(i, j, timeFactors, max(i), apertureInfo.totalNumOfShapes);
    timeBNOptAngles = timeFacMatrix * timeDAOBorderAngles;

    % prep
    leftLeafDiff    = zeros(apertureInfo.arc.numLeafSpeedConstraint * apertureInfo.beam(1).numOfActiveLeafPairs, 1);
    rightLeafDiff   = zeros(apertureInfo.arc.numLeafSpeedConstraint * apertureInfo.beam(1).numOfActiveLeafPairs, 1);
    tVec            = zeros(apertureInfo.arc.numLeafSpeedConstraint * apertureInfo.beam(1).numOfActiveLeafPairs, 1);
    maxLeafSpeed    = zeros(1, max(i));

    offset      = 0;
    shapeInd    = 1;

    for i = 1:numel(apertureInfo.beam)
        % loop over beams
        n = apertureInfo.beam(i).numOfActiveLeafPairs;

        if ~isempty(apertureInfo.arc.beam(i).leafConstMask)

            % get vector indices
            if apertureInfo.arc.beam(i).isDAOBeam
                % if it's a DAO beam, use own vector offset
                vectorIx_LI = apertureInfo.beam(i).shape(1).vectorOffset(1) + ((1:n) - 1);
                vectorIx_LF = apertureInfo.beam(i).shape(1).vectorOffset(2) + ((1:n) - 1);
            else
                % otherwise, use vector offset of previous and next
                % beams
                vectorIx_LI = apertureInfo.beam(apertureInfo.arc.beam(i).lastDAOBeamIx).shape(1).vectorOffset(2) + ((1:n) - 1);
                vectorIx_LF = apertureInfo.beam(apertureInfo.arc.beam(i).nextDAOBeamIx).shape(1).vectorOffset(1) + ((1:n) - 1);
            end
            vectorIx_RI = vectorIx_LI + apertureInfo.totalNumOfLeafPairs;
            vectorIx_RF = vectorIx_LF + apertureInfo.totalNumOfLeafPairs;

            % extract leaf positions, time
            leftLeafPosInitial   = apertureInfoVec(vectorIx_LI);
            rightLeafPosInitial  = apertureInfoVec(vectorIx_RI);
            leftLeafPosFinal   = apertureInfoVec(vectorIx_LF);
            rightLeafPosFinal  = apertureInfoVec(vectorIx_RF);
            t               = timeBNOptAngles(shapeInd);

            % determine indices
            indInDiffVec = offset + (1:n);

            % insert differences, time
            leftLeafDiff(indInDiffVec)  = abs(leftLeafPosFinal - leftLeafPosInitial);
            rightLeafDiff(indInDiffVec) = abs(rightLeafPosFinal - rightLeafPosInitial);
            tVec(indInDiffVec)          = t;

            % get max speed
            leftLeafSpeed = abs(leftLeafPosFinal - leftLeafPosInitial) ./ t;
            rightLeafSpeed = abs(rightLeafPosFinal - rightLeafPosInitial) ./ t;
            maxLeafSpeed_temp = max([leftLeafSpeed; rightLeafSpeed]);

            % update max speed
            if maxLeafSpeed_temp > maxLeafSpeed(shapeInd)
                maxLeafSpeed(shapeInd) = maxLeafSpeed_temp;
            end

            % update offset
            offset = offset + n;

            % increment shapeInd only for beams which have transition
            % defined
            shapeInd = shapeInd + 1;
        end
    end
else
    % value of constraints for leaves
    leftLeafPos  = apertureInfoVec((1:(apertureInfo.totalNumOfLeafPairs)) + apertureInfo.totalNumOfShapes);
    rightLeafIx  = (1 + apertureInfo.totalNumOfLeafPairs + apertureInfo.totalNumOfShapes): ...
        (apertureInfo.totalNumOfShapes + apertureInfo.totalNumOfLeafPairs * 2);
    rightLeafPos = apertureInfoVec(rightLeafIx);

    % Using the static fluence calculation, we have the leaf positions in
    % the vector be the leaf positions at the centre of the Dij arcs (for optimized angles only).
    % Therefore we must use the times between the centres of the Dij arcs (for optimized angles only).
    i = sort(repmat(1:(apertureInfo.totalNumOfShapes - 1), 1, 2));
    j = sort(repmat(1:apertureInfo.totalNumOfShapes, 1, 2));
    j(1) = [];
    j(end) = [];

    timeFactors = [apertureInfo.arc.beam.timeFactors]';
    timeFactors(1) = [];
    timeFactors(end) = [];

    timeFacMatrix = sparse(i, j, timeFactors, apertureInfo.totalNumOfShapes - 1, apertureInfo.totalNumOfShapes);
    timeBNOptAngles = timeFacMatrix * timeDAOBorderAngles;

    leftLeafSpeed = abs(diff(reshape(leftLeafPos, apertureInfo.beam(1).numOfActiveLeafPairs, []), 1, 2)) ./ ...
        repmat(timeBNOptAngles', apertureInfo.beam(1).numOfActiveLeafPairs, 1);
    rightLeafSpeed = abs(diff(reshape(rightLeafPos, apertureInfo.beam(1).numOfActiveLeafPairs, []), 1, 2)) ./ ...
        repmat(timeBNOptAngles', apertureInfo.beam(1).numOfActiveLeafPairs, 1);

    % values of max leaf speeds
    leftMaxLeafSpeed = max(leftLeafSpeed, [], 1);
    rightMaxLeafSpeed = max(rightLeafSpeed, [], 1);
    maxLeafSpeed = max([leftMaxLeafSpeed; rightMaxLeafSpeed], [], 1);
end

% enter into apertureInfo
l = 1;
maxMaxLeafSpeed = 0;
for i = 1:size(apertureInfo.beam, 2)
    if apertureInfo.arc.beam(i).isDAOBeam
        if apertureInfo.continuousAperture
            % for dynamic, we take the max leaf speed to be the actual leaf
            % speed
            ind = apertureInfo.arc.beam(i).timeFactorIx(apertureInfo.arc.beam(i).timeFactors ~= 0);

            apertureInfo.beam(i).maxLeafSpeed = max(maxLeafSpeed(ind));
            if apertureInfo.beam(i).maxLeafSpeed >= maxMaxLeafSpeed
                maxMaxLeafSpeed = apertureInfo.beam(i).maxLeafSpeed;
            end
        else
            % for static, we take the max leaf speed to be the max leaf
            % of two speeds, one being the speed in the first half-arc, the
            % second being the speed in the second half-arc (these will be
            % different in general)

            if l == 1
                apertureInfo.beam(i).maxLeafSpeed = maxLeafSpeed(l);
            elseif l == apertureInfo.totalNumOfShapes
                apertureInfo.beam(i).maxLeafSpeed = maxLeafSpeed(l - 1);
            else
                apertureInfo.beam(i).maxLeafSpeed = max(maxLeafSpeed(l - 1), maxLeafSpeed(l));
            end

            if l < apertureInfo.totalNumOfShapes && maxLeafSpeed(l) >= maxMaxLeafSpeed
                maxMaxLeafSpeed = maxLeafSpeed(l);
            end
        end

        if l < apertureInfo.totalNumOfShapes && maxLeafSpeed(l) >= maxMaxLeafSpeed
            maxMaxLeafSpeed = maxLeafSpeed(l);
        end

        l = l + 1;
    end
end

apertureInfo.maxLeafSpeed = maxMaxLeafSpeed;
