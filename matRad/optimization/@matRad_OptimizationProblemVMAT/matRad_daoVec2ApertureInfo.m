function updatedInfo = matRad_daoVec2ApertureInfo(apertureInfo, apertureInfoVect)
% matRad function to translate the vector representation of the aperture
% shape and weight into an aperture info struct. At the same time, the
% updated bixel weight vector w is computed and the analytic Jacobian
% between bixel weights and aperture vector entries is assembled
%
% call
%   updatedInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVect)
%
% input
%   apertureInfo:     aperture shape info struct
%   apertureInfoVect: aperture weights and shapes parameterized as vector
%
% output
%   updatedInfo: updated aperture shape info struct according to apertureInfoVect
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function to update the apertureInfo struct after the each iteration of the
% optimization

% This is the VMAT version: one shape per beam, arc interpolation of
% non-DAO beams, gantry time entries at the end of the vector, and
% assembly of the analytic bixel-aperture Jacobian (bixelJApVec). Static
% (step-and-shoot) plans must use matRad_OptimizationProblemDAO's version.
if ~apertureInfo.runVMAT
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError(['Static apertureInfo passed to the VMAT vector conversion! ' ...
                          'Use matRad_OptimizationProblemDAO.matRad_daoVec2ApertureInfo (or instance dispatch) instead.']);
end

% initializing variables
updatedInfo = apertureInfo;

updatedInfo.apertureVector = apertureInfoVect;

% options for bixel and Jacobian calculation
mlcOptions.bixelWidth = apertureInfo.bixelWidth;
calcOptions.continuousAperture = updatedInfo.continuousAperture;
calcOptions.saveJacobian = true;
vectorIndices.totalNumOfShapes = apertureInfo.totalNumOfShapes;

w = zeros(apertureInfo.totalNumOfBixels, 1);

% first index of the gantry time entries in the aperture vector
tIxOffset = updatedInfo.totalNumOfShapes + updatedInfo.totalNumOfLeafPairs * 2;

if ~all([updatedInfo.arc.beam.isDAOBeam])
    % pre-pass: update weights/times of all DAO beams first, since the
    % interpolated beams below reference their last/next DAO beam - which
    % may come later in the main loop
    j = 1;
    for i = 1:numel(updatedInfo.beam)
        if updatedInfo.arc.beam(i).isDAOBeam
            % update the shape weight
            % rescale the weight from the vector using the previous
            % iteration scaling factor
            updatedInfo.beam(i).shape(j).weight = apertureInfoVect(updatedInfo.beam(i).shape(j).weightOffset) ./ ...
                updatedInfo.beam(i).shape(j).jacobiScale;

            updatedInfo.beam(i).shape(j).MU = updatedInfo.beam(i).shape(j).weight * updatedInfo.weightToMU;
            updatedInfo.beam(i).time = apertureInfoVect(tIxOffset + updatedInfo.arc.beam(i).DAOBeamNumber) * ...
                updatedInfo.arc.beam(i).timeFactorCurrent;
            updatedInfo.beam(i).gantryRot = updatedInfo.arc.beam(i).doseAngleBordersDiff / updatedInfo.beam(i).time;
            updatedInfo.beam(i).shape(j).MURate = updatedInfo.beam(i).shape(j).MU ./ updatedInfo.beam(i).time;
        end
    end
end

% Jacobian matrix to be used in the DAO gradient function
% this tells us the gradient of a particular bixel with respect to an
% element in the apertureVector (aperture weight or leaf position)
% store as a vector for now, convert to sparse matrix later

optBixelFactor = 7;
% For optimized beams: 7 = (1 from weights) + (3 from left leaf positions (I, M, and F)) + (3 from
% right leaf positions (I, M, and F))

intBixelFactor = 2 * optBixelFactor + 2;
% For interpolated beams: multiply this number times 2 (influenced by the
% one before and the one after), then add 2 (influenced by the time of the
% times before and after)

% for the time (probability) gradients
optBixelFactor = optBixelFactor + apertureInfo.totalNumOfShapes;
intBixelFactor = intBixelFactor + apertureInfo.totalNumOfShapes;

bixelJApVec_sz = (updatedInfo.totalNumOfOptBixels * optBixelFactor + ...
                  (updatedInfo.totalNumOfBixels - updatedInfo.totalNumOfOptBixels) * intBixelFactor) * 2;

bixelJApVecVec = zeros(1, bixelJApVec_sz);

% vector indices
bixelJApVecI = nan(1, bixelJApVec_sz);
% bixel indices
bixelJApVecJ = zeros(1, bixelJApVec_sz);
% offset
bixelJApVecOffset = 0;

%% update the shapeMaps
% here the new colimator positions are used to create new shapeMaps that
% now include decimal values instead of binary

% loop over all beams
for i = 1:numel(updatedInfo.beam)

    % pre compute left and right bixel edges
    edgesLeft = updatedInfo.beam(i).posOfCornerBixel(1) + ...
        ((1:size(apertureInfo.beam(i).bixelIndMap, 2)) - 1 - 1 / 2) * updatedInfo.bixelWidth;
    edgesRight = updatedInfo.beam(i).posOfCornerBixel(1) + ...
        ((1:size(apertureInfo.beam(i).bixelIndMap, 2)) - 1 + 1 / 2) * updatedInfo.bixelWidth;

    % get dimensions of 2d matrices that store shape/bixel information
    n = apertureInfo.beam(i).numOfActiveLeafPairs;

    % in VMAT there is always exactly one shape per beam
    numOfShapes = 1;
    calcOptions.isDAOBeam = updatedInfo.arc.beam(i).isDAOBeam;

    mlcOptions.limLeft = apertureInfo.beam(i).limLeft;
    mlcOptions.limRight = apertureInfo.beam(i).limRight;
    mlcOptions.edgesLeft = edgesLeft;
    mlcOptions.edgesRight = edgesRight;
    mlcOptions.centres = (edgesLeft + edgesRight) / 2;
    mlcOptions.widths = edgesRight - edgesLeft;
    mlcOptions.n = n;
    mlcOptions.numBix = size(apertureInfo.beam(i).bixelIndMap, 2);
    mlcOptions.bixelIndMap = apertureInfo.beam(i).bixelIndMap;

    for j = 1:numOfShapes

        if updatedInfo.arc.beam(i).isDAOBeam
            % this is a DAO beam

            % update the shape weight
            updatedInfo.beam(i).shape(j).weight = apertureInfoVect(updatedInfo.beam(i).shape(j).weightOffset) ./ ...
                updatedInfo.beam(i).shape(j).jacobiScale;

            updatedInfo.beam(i).shape(j).MU = updatedInfo.beam(i).shape(j).weight * updatedInfo.weightToMU;
            updatedInfo.beam(i).time = apertureInfoVect(tIxOffset + updatedInfo.arc.beam(i).DAOBeamNumber) * ...
                updatedInfo.arc.beam(i).timeFactorCurrent;
            updatedInfo.beam(i).gantryRot = updatedInfo.arc.beam(i).doseAngleBordersDiff / updatedInfo.beam(i).time;
            updatedInfo.beam(i).shape(j).MURate = updatedInfo.beam(i).shape(j).MU ./ updatedInfo.beam(i).time;

            if ~updatedInfo.continuousAperture
                % extract left and right leaf positions from shape vector
                vectorIx_L = updatedInfo.beam(i).shape(j).vectorOffset + ((1:n) - 1);
                vectorIx_R = vectorIx_L + apertureInfo.totalNumOfLeafPairs;
                leftLeafPos  = apertureInfoVect(vectorIx_L);
                rightLeafPos = apertureInfoVect(vectorIx_R);

                % update information in shape structure
                updatedInfo.beam(i).shape(j).leftLeafPos  = leftLeafPos;
                updatedInfo.beam(i).shape(j).leftLeafPosInitial = leftLeafPos;
                updatedInfo.beam(i).shape(j).leftLeafPosFinal = leftLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPos = rightLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPosInitial = rightLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPosFinal = rightLeafPos;
            else
                % extract left and right leaf positions from shape vector
                vectorIxLI = updatedInfo.beam(i).shape(j).vectorOffset(1) + ((1:n) - 1);
                vectorIxRI = vectorIxLI + apertureInfo.totalNumOfLeafPairs;
                leftLeafPosInitial = apertureInfoVect(vectorIxLI);
                rightLeafPosInitial = apertureInfoVect(vectorIxRI);

                vectorIxLF = updatedInfo.beam(i).shape(j).vectorOffset(2) + ((1:n) - 1);
                vectorIxRF = vectorIxLF + apertureInfo.totalNumOfLeafPairs;
                leftLeafPosFinal = apertureInfoVect(vectorIxLF);
                rightLeafPosFinal = apertureInfoVect(vectorIxRF);

                % update information in shape structure
                updatedInfo.beam(i).shape(j).leftLeafPosInitial  = leftLeafPosInitial;
                updatedInfo.beam(i).shape(j).rightLeafPosInitial = rightLeafPosInitial;

                updatedInfo.beam(i).shape(j).leftLeafPosFinal  = leftLeafPosFinal;
                updatedInfo.beam(i).shape(j).rightLeafPosFinal = rightLeafPosFinal;
            end

        else
            % this is an interpolated beam
            lastDAOIx = updatedInfo.arc.beam(i).lastDAOBeamIx;
            nextDAOIx = updatedInfo.arc.beam(i).nextDAOBeamIx;
            weightFracFromLastDAO = updatedInfo.arc.beam(i).weightFracFromLastDAO;

            % MURate is interpolated between MURates of optimized apertures
            updatedInfo.beam(i).gantryRot = 1 ./ (updatedInfo.arc.beam(i).timeFracFromLastDAO ./ updatedInfo.beam(lastDAOIx).gantryRot + ...
                                                  updatedInfo.arc.beam(i).timeFracFromNextDAO ./ updatedInfo.beam(nextDAOIx).gantryRot);
            updatedInfo.beam(i).time = updatedInfo.arc.beam(i).doseAngleBordersDiff ./ updatedInfo.beam(i).gantryRot;
            updatedInfo.beam(i).shape(j).MURate = weightFracFromLastDAO * updatedInfo.beam(lastDAOIx).shape(j).MURate + ...
                (1 - weightFracFromLastDAO) * updatedInfo.beam(nextDAOIx).shape(j).MURate;

            % calculate MU, weight
            updatedInfo.beam(i).shape(j).MU = updatedInfo.beam(i).shape(j).MURate .* updatedInfo.beam(i).time;
            updatedInfo.beam(i).shape(j).weight = updatedInfo.beam(i).shape(j).MU ./ updatedInfo.weightToMU;

            if ~updatedInfo.continuousAperture

                fracFromLastOpt = updatedInfo.arc.beam(i).weightFracFromLastDAO;
                fracFromLastOptI = updatedInfo.arc.beam(i).weightFracFromLastDAO * ones(n, 1);
                fracFromLastOptF = updatedInfo.arc.beam(i).weightFracFromLastDAO * ones(n, 1);
                fracFromNextOptI = (1 - updatedInfo.arc.beam(i).weightFracFromLastDAO) * ones(n, 1);
                fracFromNextOptF = (1 - updatedInfo.arc.beam(i).weightFracFromLastDAO) * ones(n, 1);

                % obtain leaf positions at last DAO beam
                vectorIxLFLast = updatedInfo.beam(lastDAOIx).shape(j).vectorOffset + ((1:n) - 1);
                vectorIxRFLast = vectorIxLFLast + apertureInfo.totalNumOfLeafPairs;
                leftLeafPos_last = apertureInfoVect(vectorIxLFLast);
                rightLeafPos_last = apertureInfoVect(vectorIxRFLast);

                % obtain leaf positions at next DAO beam
                vectorIxLINext = updatedInfo.beam(nextDAOIx).shape(j).vectorOffset + ((1:n) - 1);
                vectorIxRINext = vectorIxLINext + apertureInfo.totalNumOfLeafPairs;
                leftLeafPos_next = apertureInfoVect(vectorIxLINext);
                rightLeafPos_next = apertureInfoVect(vectorIxRINext);

                % interpolate leaf positions
                leftLeafPos = weightFracFromLastDAO * leftLeafPos_last + (1 - weightFracFromLastDAO) * leftLeafPos_next;
                rightLeafPos = weightFracFromLastDAO * rightLeafPos_last + (1 - weightFracFromLastDAO) * rightLeafPos_next;

                % update information in shape structure
                updatedInfo.beam(i).shape(j).leftLeafPos  = leftLeafPos;
                updatedInfo.beam(i).shape(j).leftLeafPosInitial = leftLeafPos;
                updatedInfo.beam(i).shape(j).leftLeafPosFinal = leftLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPos = rightLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPosInitial = rightLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPosFinal = rightLeafPos;
            else

                fracFromLastOpt = updatedInfo.arc.beam(i).weightFracFromLastDAO;
                fracFromLastOptI = updatedInfo.arc.beam(i).weightFracFromLastDAOInitial * ones(n, 1);
                fracFromLastOptF = updatedInfo.arc.beam(i).weightFracFromLastDAOFinal * ones(n, 1);
                fracFromNextOptI = updatedInfo.arc.beam(i).weightFracFromNextDAOInitial * ones(n, 1);
                fracFromNextOptF = updatedInfo.arc.beam(i).weightFracFromNextDAOFinal * ones(n, 1);

                % obtain leaf positions at last DAO beam
                vectorIxLFLast = updatedInfo.beam(updatedInfo.arc.beam(i).lastDAOBeamIx).shape(j).vectorOffset(2) + ((1:n) - 1);
                vectorIxRFLast = vectorIxLFLast + apertureInfo.totalNumOfLeafPairs;
                leftLeafPos_F_last = apertureInfoVect(vectorIxLFLast);
                rightLeafPos_F_last = apertureInfoVect(vectorIxRFLast);

                % obtain leaf positions at next DAO beam
                vectorIxLINext = updatedInfo.beam(updatedInfo.arc.beam(i).nextDAOBeamIx).shape(j).vectorOffset(1) + ((1:n) - 1);
                vectorIxRINext = vectorIxLINext + apertureInfo.totalNumOfLeafPairs;
                leftLeafPos_I_next = apertureInfoVect(vectorIxLINext);
                rightLeafPos_I_next = apertureInfoVect(vectorIxRINext);

                % interpolate leaf positions
                updatedInfo.beam(i).shape(j).leftLeafPosInitial = fracFromLastOptI .* leftLeafPos_F_last + fracFromNextOptI .* leftLeafPos_I_next;
                updatedInfo.beam(i).shape(j).rightLeafPosInitial = fracFromLastOptI .* rightLeafPos_F_last + fracFromNextOptI .* rightLeafPos_I_next;

                updatedInfo.beam(i).shape(j).leftLeafPosFinal = fracFromLastOptF .* leftLeafPos_F_last + fracFromNextOptF .* leftLeafPos_I_next;
                updatedInfo.beam(i).shape(j).rightLeafPosFinal = fracFromLastOptF .* rightLeafPos_F_last + fracFromNextOptF .* rightLeafPos_I_next;
            end
        end

    end

    for j = 1:numOfShapes

        % shapeMap
        shapeMap = zeros(size(updatedInfo.beam(i).bixelIndMap));
        % sumGradSq
        sumGradSq = 0;

        % insert variables
        vectorIndices.tIx_Vec       = (apertureInfo.totalNumOfShapes + apertureInfo.totalNumOfLeafPairs * 2) + (1:apertureInfo.totalNumOfShapes);

        variables.weight            = updatedInfo.beam(i).shape(j).weight;
        variables.leftLeafPosInitial     = updatedInfo.beam(i).shape(j).leftLeafPosInitial;
        variables.leftLeafPosFinal     = updatedInfo.beam(i).shape(j).leftLeafPosFinal;
        variables.rightLeafPosInitial    = updatedInfo.beam(i).shape(j).rightLeafPosInitial;
        variables.rightLeafPosFinal    = updatedInfo.beam(i).shape(j).rightLeafPosFinal;

        if updatedInfo.arc.beam(i).isDAOBeam

            variables.jacobiScale = updatedInfo.beam(i).shape(1).jacobiScale;

            vectorIndices.DAOBeamNumber      = updatedInfo.arc.beam(i).DAOBeamNumber;
            if updatedInfo.continuousAperture
                vectorIndices.vectorIxLI   = updatedInfo.beam(i).shape(j).vectorOffset(1) + ((1:n) - 1);
                vectorIndices.vectorIxLF   = updatedInfo.beam(i).shape(j).vectorOffset(2) + ((1:n) - 1);
                vectorIndices.vectorIxRI   = vectorIndices.vectorIxLI + apertureInfo.totalNumOfLeafPairs;
                vectorIndices.vectorIxRF   = vectorIndices.vectorIxLF + apertureInfo.totalNumOfLeafPairs;
            else
                vectorIndices.vectorIxLI   = updatedInfo.beam(i).shape(j).vectorOffset + ((1:n) - 1);
                vectorIndices.vectorIxLF   = updatedInfo.beam(i).shape(j).vectorOffset + ((1:n) - 1);
                vectorIndices.vectorIxRI   = vectorIndices.vectorIxLI + apertureInfo.totalNumOfLeafPairs;
                vectorIndices.vectorIxRF   = vectorIndices.vectorIxLF + apertureInfo.totalNumOfLeafPairs;
            end
        else

            variables.weightLast = updatedInfo.beam(updatedInfo.arc.beam(i).lastDAOBeamIx).shape(j).weight;
            variables.weightNext = updatedInfo.beam(updatedInfo.arc.beam(i).nextDAOBeamIx).shape(j).weight;

            variables.jacobiScaleLast    = updatedInfo.beam(updatedInfo.arc.beam(i).lastDAOBeamIx).shape(1).jacobiScale;
            variables.jacobiScaleNext    = updatedInfo.beam(updatedInfo.arc.beam(i).nextDAOBeamIx).shape(1).jacobiScale;

            variables.timeLast = updatedInfo.beam(updatedInfo.arc.beam(i).lastDAOBeamIx).time;
            variables.timeNext = updatedInfo.beam(updatedInfo.arc.beam(i).nextDAOBeamIx).time;
            variables.time      = updatedInfo.beam(i).time;

            variables.fracFromLastOptI  = fracFromLastOptI;
            variables.fracFromLastOptF  = fracFromLastOptF;
            variables.fracFromNextOptI  = fracFromNextOptI;
            variables.fracFromNextOptF  = fracFromNextOptF;
            variables.fracFromLastOpt   = fracFromLastOpt;

            variables.doseAngleBordersDiff      = updatedInfo.arc.beam(i).doseAngleBordersDiff;
            variables.doseAngleBordersDiffLast = updatedInfo.arc.beam(updatedInfo.arc.beam(i).lastDAOBeamIx).doseAngleBordersDiff;
            variables.doseAngleBordersDiffNext = updatedInfo.arc.beam(updatedInfo.arc.beam(i).nextDAOBeamIx).doseAngleBordersDiff;
            variables.timeFactorCurrentLast          = updatedInfo.arc.beam(updatedInfo.arc.beam(i).lastDAOBeamIx).timeFactorCurrent;
            variables.timeFactorCurrentNext          = updatedInfo.arc.beam(updatedInfo.arc.beam(i).nextDAOBeamIx).timeFactorCurrent;
            variables.weightFracFromLastDAO           = updatedInfo.arc.beam(i).weightFracFromLastDAO;
            variables.timeFracFromLastDAO       = updatedInfo.arc.beam(i).timeFracFromLastDAO;
            variables.timeFracFromNextDAO       = updatedInfo.arc.beam(i).timeFracFromNextDAO;

            vectorIndices.DAOBeamNumberLast = updatedInfo.arc.beam(updatedInfo.arc.beam(i).lastDAOBeamIx).DAOBeamNumber;
            vectorIndices.DAOBeamNumberNext = updatedInfo.arc.beam(updatedInfo.arc.beam(i).nextDAOBeamIx).DAOBeamNumber;
            vectorIndices.tIxLast      = tIxOffset + vectorIndices.DAOBeamNumberLast;
            vectorIndices.tIxNext      = tIxOffset + vectorIndices.DAOBeamNumberNext;

            if updatedInfo.continuousAperture
                vectorIndices.vectorIxLFLast  = updatedInfo.beam(updatedInfo.arc.beam(i).lastDAOBeamIx).shape(j).vectorOffset(2) + ((1:n) - 1);
                vectorIndices.vectorIxLINext  = updatedInfo.beam(updatedInfo.arc.beam(i).nextDAOBeamIx).shape(j).vectorOffset(1) + ((1:n) - 1);
                vectorIndices.vectorIxRFLast  = vectorIndices.vectorIxLFLast + apertureInfo.totalNumOfLeafPairs;
                vectorIndices.vectorIxRINext  = vectorIndices.vectorIxLINext + apertureInfo.totalNumOfLeafPairs;
            else
                vectorIndices.vectorIxLFLast  = updatedInfo.beam(updatedInfo.arc.beam(i).lastDAOBeamIx).shape(j).vectorOffset + ((1:n) - 1);
                vectorIndices.vectorIxLINext  = updatedInfo.beam(updatedInfo.arc.beam(i).nextDAOBeamIx).shape(j).vectorOffset + ((1:n) - 1);
                vectorIndices.vectorIxRFLast  = vectorIndices.vectorIxLFLast + apertureInfo.totalNumOfLeafPairs;
                vectorIndices.vectorIxRINext  = vectorIndices.vectorIxLINext + apertureInfo.totalNumOfLeafPairs;
            end
        end

        counters.bixelJApVecOffset = bixelJApVecOffset;

        % calculate bixel weight and derivative in function
        accum.w               = w;
        accum.bixelJApVec.vec = bixelJApVecVec;
        accum.bixelJApVec.i   = bixelJApVecI;
        accum.bixelJApVec.j   = bixelJApVecJ;
        accum.sumGradSq       = sumGradSq;
        accum.shapeMapW       = shapeMap;
        accum.counters        = counters;

        accum = matRad_calcBixelWeightAndGradient(calcOptions, mlcOptions, variables, vectorIndices, accum);

        w              = accum.w;
        bixelJApVecVec = accum.bixelJApVec.vec;
        bixelJApVecI   = accum.bixelJApVec.i;
        bixelJApVecJ   = accum.bixelJApVec.j;
        sumGradSq      = accum.sumGradSq;
        shapeMap       = accum.shapeMapW;
        counters       = accum.counters;

        bixelJApVecOffset = counters.bixelJApVecOffset;

        % update shapeMap
        updatedInfo.beam(i).shape(j).shapeMap = shapeMap;
        % update sumGradSq
        % FIX THIS FOR INTERPOLATED ANGLES???
        updatedInfo.beam(i).shape(j).sumGradSq = sumGradSq;

    end
end

% save bixelWeight, apertureVector
updatedInfo.bixelWeights = w;
updatedInfo.apertureVector = apertureInfoVect;

% save Jacobian between bixelWeight, apertureVector
deleteInd_i = isnan(bixelJApVecI);
deleteInd_j = bixelJApVecJ == 0;
if ~all(deleteInd_i == deleteInd_j)
    error('Jacobian deletion mismatch');
else
    bixelJApVecI(deleteInd_i) = [];
    bixelJApVecJ(deleteInd_i) = [];
    bixelJApVecVec(deleteInd_i) = [];
end
updatedInfo.bixelJApVec = sparse(bixelJApVecI, bixelJApVecJ, bixelJApVecVec, numel(apertureInfoVect), updatedInfo.totalNumOfBixels);

end
