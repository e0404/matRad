function jacob = matRad_constraintJacobian(optiProb, apertureInfoVec, dij, cst)
% matRad IPOPT callback: jacobian function for VMAT optimization
%
% call
%   jacob = matRad_constraintJacobian(optiProb,apertureInfoVec,dij,cst)
%
% input
%   apertureInfoVec: aperture info vector
%   dij:             dose influence matrix
%   cst:             matRad cst struct
%
% output
%   jacob:           jacobian of constraint function
%
% References
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

% dosimetric and leaf-gap (DAO) constraint jacobians from the superclass;
% this also updates optiProb.apertureInfo for the vector if necessary
jacob_dos_dao = matRad_constraintJacobian@matRad_OptimizationProblemDAO(optiProb, apertureInfoVec, dij, cst);
apertureInfo = optiProb.apertureInfo;

% VMAT machine constraints (leaf speed, dose rate)
% values of times spent in an arc surrounding the optimized angles (full
% arc/dose influence arc)
timeDAOBorderAngles = apertureInfoVec(((apertureInfo.totalNumOfShapes + apertureInfo.totalNumOfLeafPairs * 2) + 1):end);
timeFactorCurrent = [apertureInfo.arc.beam([apertureInfo.arc.beam.isDAOBeam]).timeFactorCurrent]';
timeDoseBorderAngles = timeDAOBorderAngles .* timeFactorCurrent;

if apertureInfo.continuousAperture
    timeFactors = [apertureInfo.arc.beam.timeFactors]';
    deleteInd = timeFactors == 0;
    timeFactors(deleteInd) = [];

    i = [apertureInfo.arc.beam.timeFactorIx]';
    i(deleteInd) = [];

    j = repelem(1:apertureInfo.totalNumOfShapes, 1, 3);
    j(deleteInd) = [];

    timeFacMatrix = sparse(i, j, timeFactors, max(i), apertureInfo.totalNumOfShapes);
    timeBNOptAngles = timeFacMatrix * timeDAOBorderAngles;

    % set up
    n = apertureInfo.beam(1).numOfActiveLeafPairs;
    indInSparseVec  = (1:n);
    indInConVec     = (1:n);
    shapeInd        = 1;

    % sparse matrix
    numElem     = n .* (apertureInfo.arc.numLeafSpeedConstraintDAO * 6 + ...
                        (apertureInfo.arc.numLeafSpeedConstraint - apertureInfo.arc.numLeafSpeedConstraintDAO) * 8);
    i_sparse    = zeros(numElem, 1);
    j_sparse    = zeros(numElem, 1);
    s_sparse    = zeros(numElem, 1);

    for i = 1:numel(apertureInfo.beam)
        % loop over beams

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
            leftLeafPos_I   = apertureInfoVec(vectorIx_LI);
            rightLeafPos_I  = apertureInfoVec(vectorIx_RI);
            leftLeafPos_F   = apertureInfoVec(vectorIx_LF);
            rightLeafPos_F  = apertureInfoVec(vectorIx_RF);
            t               = timeBNOptAngles(shapeInd);

            % calc diffs
            leftLeafDiff    = leftLeafPos_F - leftLeafPos_I;
            rightLeafDiff   = rightLeafPos_F - rightLeafPos_I;

            % calc jacobs

            % wrt initial leaf positions (left, then right)
            i_sparse(indInSparseVec)    = indInConVec;
            j_sparse(indInSparseVec)    = vectorIx_LI;
            s_sparse(indInSparseVec)    = -sign(leftLeafDiff) ./ t;
            indInSparseVec              = indInSparseVec + n;

            i_sparse(indInSparseVec)    = indInConVec + apertureInfo.arc.numLeafSpeedConstraint * apertureInfo.beam(1).numOfActiveLeafPairs;
            j_sparse(indInSparseVec)    = vectorIx_RI;
            s_sparse(indInSparseVec)    = -sign(rightLeafDiff) ./ t;
            indInSparseVec              = indInSparseVec + n;

            % wrt final leaf positions (left, then right)
            i_sparse(indInSparseVec)    = indInConVec;
            j_sparse(indInSparseVec)    = vectorIx_LF;
            s_sparse(indInSparseVec)    = sign(leftLeafDiff) ./ t;
            indInSparseVec              = indInSparseVec + n;

            i_sparse(indInSparseVec)    = indInConVec + apertureInfo.arc.numLeafSpeedConstraint * apertureInfo.beam(1).numOfActiveLeafPairs;
            j_sparse(indInSparseVec)    = vectorIx_RF;
            s_sparse(indInSparseVec)    = sign(rightLeafDiff) ./ t;
            indInSparseVec              = indInSparseVec + n;

            % wrt time (left, then right)
            % how we do this depends on if it's a DAO beam or
            % not
            if apertureInfo.arc.beam(i).isDAOBeam
                % if it is, then speeds only depend on its own
                % time
                i_sparse(indInSparseVec)    = indInConVec;
                j_sparse(indInSparseVec)    = apertureInfo.arc.beam(i).timeInd;
                s_sparse(indInSparseVec)    = -apertureInfo.arc.beam(i).timeFactors(2) .* abs(leftLeafDiff) ./ (t.^2);
                indInSparseVec              = indInSparseVec + n;

                i_sparse(indInSparseVec)    = indInConVec + apertureInfo.arc.numLeafSpeedConstraint * apertureInfo.beam(1).numOfActiveLeafPairs;
                j_sparse(indInSparseVec)    = apertureInfo.arc.beam(i).timeInd;
                s_sparse(indInSparseVec)    = -apertureInfo.arc.beam(i).timeFactors(2) .* abs(rightLeafDiff) ./ (t.^2);
                indInSparseVec              = indInSparseVec + n;

            else
                % otherwise, speed depends on time of DAO
                % before and DAO after
                lastDAOBeam = apertureInfo.arc.beam(apertureInfo.arc.beam(i).lastDAOBeamIx);
                nextDAOBeam = apertureInfo.arc.beam(apertureInfo.arc.beam(i).nextDAOBeamIx);

                % before
                i_sparse(indInSparseVec)    = indInConVec;
                j_sparse(indInSparseVec)    = lastDAOBeam.timeInd;
                s_sparse(indInSparseVec)    = -lastDAOBeam.timeFactors(3) .* abs(leftLeafDiff) ./ (t.^2);
                indInSparseVec              = indInSparseVec + n;

                i_sparse(indInSparseVec)    = indInConVec + apertureInfo.arc.numLeafSpeedConstraint * apertureInfo.beam(1).numOfActiveLeafPairs;
                j_sparse(indInSparseVec)    = lastDAOBeam.timeInd;
                s_sparse(indInSparseVec)    = -lastDAOBeam.timeFactors(3) .* abs(rightLeafDiff) ./ (t.^2);
                indInSparseVec              = indInSparseVec + n;

                % after
                i_sparse(indInSparseVec)    = indInConVec;
                j_sparse(indInSparseVec)    = nextDAOBeam.timeInd;
                s_sparse(indInSparseVec)    = -nextDAOBeam.timeFactors(1) .* abs(leftLeafDiff) ./ (t.^2);
                indInSparseVec              = indInSparseVec + n;

                i_sparse(indInSparseVec)    = indInConVec + apertureInfo.arc.numLeafSpeedConstraint * apertureInfo.beam(1).numOfActiveLeafPairs;
                j_sparse(indInSparseVec)    = nextDAOBeam.timeInd;
                s_sparse(indInSparseVec)    = -nextDAOBeam.timeFactors(1) .* abs(rightLeafDiff) ./ (t.^2);
                indInSparseVec              = indInSparseVec + n;

            end

            % update offset
            indInConVec = indInConVec + n;

            % increment shapeInd only for beams which have transition
            % defined
            shapeInd = shapeInd + 1;
        end
    end

    jacob_lfspd = sparse(i_sparse, j_sparse, s_sparse, ...
                         2 * apertureInfo.beam(1).numOfActiveLeafPairs * apertureInfo.arc.numLeafSpeedConstraint, ...
                         numel(apertureInfoVec));

else

    % shorthand aliases
    nShapes = apertureInfo.totalNumOfShapes;
    nPairs  = apertureInfo.beam(1).numOfActiveLeafPairs;
    nLP     = apertureInfo.totalNumOfLeafPairs;

    % get index values for the jacobian
    % variable index
    % value of constraints for leaves
    leftLeafPos  = apertureInfoVec((1:nLP) + nShapes);
    rightLeafPos = apertureInfoVec(1 + nLP + nShapes:nShapes + nLP * 2);

    i = sort(repmat(1:(nShapes - 1), 1, 2));
    j = sort(repmat(1:nShapes, 1, 2));
    j(1) = [];
    j(end) = [];

    timeFactors = [apertureInfo.arc.beam([apertureInfo.arc.beam.isDAOBeam]).timeFactors]';
    timeFactors(1) = [];
    timeFactors(end) = [];

    timeFacMatrix = sparse(i, j, timeFactors, nShapes - 1, nShapes);
    timeBNOptAngles = timeFacMatrix * timeDAOBorderAngles;
    nTrans = numel(timeBNOptAngles);

    currentLeftLeafInd  = (nShapes + 1):(nShapes + nPairs * nTrans);
    currentRightLeafInd = (nShapes + nLP + 1):(nShapes + nLP + nPairs * nTrans);
    nextLeftLeafInd     = (nPairs + nShapes + 1):(nPairs + nShapes + nPairs * nTrans);
    nextRightLeafInd    = (nPairs + nShapes + nLP + 1):(nPairs + nShapes + nLP + nPairs * nTrans);
    leftTimeInd  = kron(j, ones(1, nPairs)) + nShapes + nLP * 2;
    rightTimeInd = kron(j, ones(1, nPairs)) + nShapes + nLP * 2;
    % constraint index
    constraintInd = 1:2 * nPairs * nTrans;

    % jacobian of the leafspeed constraint
    i = repmat((i' - 1) * nPairs, 1, nPairs) + repmat(1:nPairs, 2 * nTrans, 1);
    i = reshape([i' i' + nPairs * nTrans], 1, []);

    i = [repmat(constraintInd, 1, 2) i];
    j = [currentLeftLeafInd currentRightLeafInd nextLeftLeafInd nextRightLeafInd leftTimeInd rightTimeInd];
    % first do jacob wrt current leaf position (left, right), then next leaf
    % position (left, right), then time (left, right)
    j_lfspd_cur = -reshape([sign(diff(reshape(leftLeafPos, nPairs, nShapes), 1, 2)) ...
                            sign(diff(reshape(rightLeafPos, nPairs, nShapes), 1, 2))] ./ ...
                           repmat(timeBNOptAngles', nPairs, 2), 2 * nPairs * nTrans, 1);

    j_lfspd_nxt = reshape([sign(diff(reshape(leftLeafPos, nPairs, nShapes), 1, 2)) ...
                           sign(diff(reshape(rightLeafPos, nPairs, nShapes), 1, 2))] ./ ...
                          repmat(timeBNOptAngles', nPairs, 2), 2 * nPairs * nTrans, 1);

    j_lfspd_t = -reshape([kron(abs(diff(reshape(leftLeafPos, nPairs, nShapes), 1, 2)), ones(1, 2)) .* repmat(timeFactors', nPairs, 1) ...
                          kron(abs(diff(reshape(rightLeafPos, nPairs, nShapes), 1, 2)), ones(1, 2)) .* repmat(timeFactors', nPairs, 1)] ./ ...
                         repmat(kron((timeBNOptAngles.^2)', ones(1, 2)), nPairs, 2), [], 1);

    s = [j_lfspd_cur; j_lfspd_nxt; j_lfspd_t];

    jacob_lfspd = sparse(i, j, s, 2 * nPairs * (nShapes - 1), numel(apertureInfoVec), numel(s));
end

% jacobian of the doserate constraint
% values of doserate (MU/sec) between optimized gantry angles
weights = apertureInfoVec(1:(apertureInfo.totalNumOfShapes)) ./ apertureInfo.jacobiScale;

timeIxOffset = apertureInfo.totalNumOfShapes + apertureInfo.totalNumOfLeafPairs * 2;
i = repmat(1:(apertureInfo.totalNumOfShapes), 1, 2);
j = [1:(apertureInfo.totalNumOfShapes) ...
     (timeIxOffset + 1):(timeIxOffset + apertureInfo.totalNumOfShapes)];
% first do jacob wrt weights, then wrt times

s = [apertureInfo.weightToMU ./ (timeDoseBorderAngles .* apertureInfo.jacobiScale); ...
     -apertureInfo.weightToMU .* weights .* timeFactorCurrent ./ (timeDoseBorderAngles.^2)];

jacob_dosrt = sparse(i, j, s, apertureInfo.totalNumOfShapes, numel(apertureInfoVec), 2 * apertureInfo.totalNumOfShapes);

% concatenate
jacob = [jacob_dos_dao; jacob_lfspd; jacob_dosrt];
