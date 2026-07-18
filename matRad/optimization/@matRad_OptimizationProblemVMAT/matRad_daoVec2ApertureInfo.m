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

if ~all([updatedInfo.propVMAT.beam.DAOBeam])
    % pre-pass: update weights/times of all DAO beams first, since the
    % interpolated beams below reference their last/next DAO beam - which
    % may come later in the main loop
    j = 1;
    for i = 1:numel(updatedInfo.beam)
        if updatedInfo.propVMAT.beam(i).DAOBeam
            % update the shape weight
            % rescale the weight from the vector using the previous
            % iteration scaling factor
            updatedInfo.beam(i).shape(j).weight = apertureInfoVect(updatedInfo.beam(i).shape(j).weightOffset) ./ ...
                updatedInfo.beam(i).shape(j).jacobiScale;

            updatedInfo.beam(i).shape(j).MU = updatedInfo.beam(i).shape(j).weight * updatedInfo.weightToMU;
            updatedInfo.beam(i).time = apertureInfoVect(tIxOffset + updatedInfo.propVMAT.beam(i).DAOIndex) * ...
                updatedInfo.propVMAT.beam(i).timeFacCurr;
            updatedInfo.beam(i).gantryRot = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff / updatedInfo.beam(i).time;
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

bixelJApVec_vec = zeros(1, bixelJApVec_sz);

% vector indices
bixelJApVec_i = nan(1, bixelJApVec_sz);
% bixel indices
bixelJApVec_j = zeros(1, bixelJApVec_sz);
% offset
bixelJApVec_offset = 0;

%% update the shapeMaps
% here the new colimator positions are used to create new shapeMaps that
% now include decimal values instead of binary

% loop over all beams
for i = 1:numel(updatedInfo.beam)

    % pre compute left and right bixel edges
    edges_l = updatedInfo.beam(i).posOfCornerBixel(1) + ...
        ((1:size(apertureInfo.beam(i).bixelIndMap, 2)) - 1 - 1 / 2) * updatedInfo.bixelWidth;
    edges_r = updatedInfo.beam(i).posOfCornerBixel(1) + ...
        ((1:size(apertureInfo.beam(i).bixelIndMap, 2)) - 1 + 1 / 2) * updatedInfo.bixelWidth;

    % get dimensions of 2d matrices that store shape/bixel information
    n = apertureInfo.beam(i).numOfActiveLeafPairs;

    % in VMAT there is always exactly one shape per beam
    numOfShapes = 1;
    calcOptions.DAOBeam = updatedInfo.propVMAT.beam(i).DAOBeam;

    mlcOptions.lim_l = apertureInfo.beam(i).lim_l;
    mlcOptions.lim_r = apertureInfo.beam(i).lim_r;
    mlcOptions.edges_l = edges_l;
    mlcOptions.edges_r = edges_r;
    mlcOptions.centres = (edges_l + edges_r) / 2;
    mlcOptions.widths = edges_r - edges_l;
    mlcOptions.n = n;
    mlcOptions.numBix = size(apertureInfo.beam(i).bixelIndMap, 2);
    mlcOptions.bixelIndMap = apertureInfo.beam(i).bixelIndMap;

    for j = 1:numOfShapes

        if updatedInfo.propVMAT.beam(i).DAOBeam
            % this is a DAO beam

            % update the shape weight
            updatedInfo.beam(i).shape(j).weight = apertureInfoVect(updatedInfo.beam(i).shape(j).weightOffset) ./ ...
                updatedInfo.beam(i).shape(j).jacobiScale;

            updatedInfo.beam(i).shape(j).MU = updatedInfo.beam(i).shape(j).weight * updatedInfo.weightToMU;
            updatedInfo.beam(i).time = apertureInfoVect(tIxOffset + updatedInfo.propVMAT.beam(i).DAOIndex) * ...
                updatedInfo.propVMAT.beam(i).timeFacCurr;
            updatedInfo.beam(i).gantryRot = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff / updatedInfo.beam(i).time;
            updatedInfo.beam(i).shape(j).MURate = updatedInfo.beam(i).shape(j).MU ./ updatedInfo.beam(i).time;

            if ~updatedInfo.continuousAperture
                % extract left and right leaf positions from shape vector
                vectorIx_L = updatedInfo.beam(i).shape(j).vectorOffset + ((1:n) - 1);
                vectorIx_R = vectorIx_L + apertureInfo.totalNumOfLeafPairs;
                leftLeafPos  = apertureInfoVect(vectorIx_L);
                rightLeafPos = apertureInfoVect(vectorIx_R);

                % update information in shape structure
                updatedInfo.beam(i).shape(j).leftLeafPos  = leftLeafPos;
                updatedInfo.beam(i).shape(j).leftLeafPos_I = leftLeafPos;
                updatedInfo.beam(i).shape(j).leftLeafPos_F = leftLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPos = rightLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPos_I = rightLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPos_F = rightLeafPos;
            else
                % extract left and right leaf positions from shape vector
                vectorIx_LI = updatedInfo.beam(i).shape(j).vectorOffset(1) + ((1:n) - 1);
                vectorIx_RI = vectorIx_LI + apertureInfo.totalNumOfLeafPairs;
                leftLeafPos_I = apertureInfoVect(vectorIx_LI);
                rightLeafPos_I = apertureInfoVect(vectorIx_RI);

                vectorIx_LF = updatedInfo.beam(i).shape(j).vectorOffset(2) + ((1:n) - 1);
                vectorIx_RF = vectorIx_LF + apertureInfo.totalNumOfLeafPairs;
                leftLeafPos_F = apertureInfoVect(vectorIx_LF);
                rightLeafPos_F = apertureInfoVect(vectorIx_RF);

                % update information in shape structure
                updatedInfo.beam(i).shape(j).leftLeafPos_I  = leftLeafPos_I;
                updatedInfo.beam(i).shape(j).rightLeafPos_I = rightLeafPos_I;

                updatedInfo.beam(i).shape(j).leftLeafPos_F  = leftLeafPos_F;
                updatedInfo.beam(i).shape(j).rightLeafPos_F = rightLeafPos_F;
            end

        else
            % this is an interpolated beam
            lastDAOIx = updatedInfo.propVMAT.beam(i).lastDAOIndex;
            nextDAOIx = updatedInfo.propVMAT.beam(i).nextDAOIndex;
            fracFromLastDAO = updatedInfo.propVMAT.beam(i).fracFromLastDAO;

            % MURate is interpolated between MURates of optimized apertures
            updatedInfo.beam(i).gantryRot = 1 ./ (updatedInfo.propVMAT.beam(i).timeFracFromLastDAO ./ updatedInfo.beam(lastDAOIx).gantryRot + ...
                                                  updatedInfo.propVMAT.beam(i).timeFracFromNextDAO ./ updatedInfo.beam(nextDAOIx).gantryRot);
            updatedInfo.beam(i).time = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff ./ updatedInfo.beam(i).gantryRot;
            updatedInfo.beam(i).shape(j).MURate = fracFromLastDAO * updatedInfo.beam(lastDAOIx).shape(j).MURate + ...
                (1 - fracFromLastDAO) * updatedInfo.beam(nextDAOIx).shape(j).MURate;

            % calculate MU, weight
            updatedInfo.beam(i).shape(j).MU = updatedInfo.beam(i).shape(j).MURate .* updatedInfo.beam(i).time;
            updatedInfo.beam(i).shape(j).weight = updatedInfo.beam(i).shape(j).MU ./ updatedInfo.weightToMU;

            if ~updatedInfo.continuousAperture

                fracFromLastOpt = updatedInfo.propVMAT.beam(i).fracFromLastDAO;
                fracFromLastOptI = updatedInfo.propVMAT.beam(i).fracFromLastDAO * ones(n, 1);
                fracFromLastOptF = updatedInfo.propVMAT.beam(i).fracFromLastDAO * ones(n, 1);
                fracFromNextOptI = (1 - updatedInfo.propVMAT.beam(i).fracFromLastDAO) * ones(n, 1);
                fracFromNextOptF = (1 - updatedInfo.propVMAT.beam(i).fracFromLastDAO) * ones(n, 1);

                % obtain leaf positions at last DAO beam
                vectorIx_LF_last = updatedInfo.beam(lastDAOIx).shape(j).vectorOffset + ((1:n) - 1);
                vectorIx_RF_last = vectorIx_LF_last + apertureInfo.totalNumOfLeafPairs;
                leftLeafPos_last = apertureInfoVect(vectorIx_LF_last);
                rightLeafPos_last = apertureInfoVect(vectorIx_RF_last);

                % obtain leaf positions at next DAO beam
                vectorIx_LI_next = updatedInfo.beam(nextDAOIx).shape(j).vectorOffset + ((1:n) - 1);
                vectorIx_RI_next = vectorIx_LI_next + apertureInfo.totalNumOfLeafPairs;
                leftLeafPos_next = apertureInfoVect(vectorIx_LI_next);
                rightLeafPos_next = apertureInfoVect(vectorIx_RI_next);

                % interpolate leaf positions
                leftLeafPos = fracFromLastDAO * leftLeafPos_last + (1 - fracFromLastDAO) * leftLeafPos_next;
                rightLeafPos = fracFromLastDAO * rightLeafPos_last + (1 - fracFromLastDAO) * rightLeafPos_next;

                % update information in shape structure
                updatedInfo.beam(i).shape(j).leftLeafPos  = leftLeafPos;
                updatedInfo.beam(i).shape(j).leftLeafPos_I = leftLeafPos;
                updatedInfo.beam(i).shape(j).leftLeafPos_F = leftLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPos = rightLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPos_I = rightLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPos_F = rightLeafPos;
            else

                fracFromLastOpt = updatedInfo.propVMAT.beam(i).fracFromLastDAO;
                fracFromLastOptI = updatedInfo.propVMAT.beam(i).fracFromLastDAO_I * ones(n, 1);
                fracFromLastOptF = updatedInfo.propVMAT.beam(i).fracFromLastDAO_F * ones(n, 1);
                fracFromNextOptI = updatedInfo.propVMAT.beam(i).fracFromNextDAO_I * ones(n, 1);
                fracFromNextOptF = updatedInfo.propVMAT.beam(i).fracFromNextDAO_F * ones(n, 1);

                % obtain leaf positions at last DAO beam
                vectorIx_LF_last = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape(j).vectorOffset(2) + ((1:n) - 1);
                vectorIx_RF_last = vectorIx_LF_last + apertureInfo.totalNumOfLeafPairs;
                leftLeafPos_F_last = apertureInfoVect(vectorIx_LF_last);
                rightLeafPos_F_last = apertureInfoVect(vectorIx_RF_last);

                % obtain leaf positions at next DAO beam
                vectorIx_LI_next = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape(j).vectorOffset(1) + ((1:n) - 1);
                vectorIx_RI_next = vectorIx_LI_next + apertureInfo.totalNumOfLeafPairs;
                leftLeafPos_I_next = apertureInfoVect(vectorIx_LI_next);
                rightLeafPos_I_next = apertureInfoVect(vectorIx_RI_next);

                % interpolate leaf positions
                updatedInfo.beam(i).shape(j).leftLeafPos_I = fracFromLastOptI .* leftLeafPos_F_last + fracFromNextOptI .* leftLeafPos_I_next;
                updatedInfo.beam(i).shape(j).rightLeafPos_I = fracFromLastOptI .* rightLeafPos_F_last + fracFromNextOptI .* rightLeafPos_I_next;

                updatedInfo.beam(i).shape(j).leftLeafPos_F = fracFromLastOptF .* leftLeafPos_F_last + fracFromNextOptF .* leftLeafPos_I_next;
                updatedInfo.beam(i).shape(j).rightLeafPos_F = fracFromLastOptF .* rightLeafPos_F_last + fracFromNextOptF .* rightLeafPos_I_next;
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
        variables.leftLeafPos_I     = updatedInfo.beam(i).shape(j).leftLeafPos_I;
        variables.leftLeafPos_F     = updatedInfo.beam(i).shape(j).leftLeafPos_F;
        variables.rightLeafPos_I    = updatedInfo.beam(i).shape(j).rightLeafPos_I;
        variables.rightLeafPos_F    = updatedInfo.beam(i).shape(j).rightLeafPos_F;

        if updatedInfo.propVMAT.beam(i).DAOBeam

            variables.jacobiScale = updatedInfo.beam(i).shape(1).jacobiScale;

            vectorIndices.DAOindex      = updatedInfo.propVMAT.beam(i).DAOIndex;
            if updatedInfo.continuousAperture
                vectorIndices.vectorIx_LI   = updatedInfo.beam(i).shape(j).vectorOffset(1) + ((1:n) - 1);
                vectorIndices.vectorIx_LF   = updatedInfo.beam(i).shape(j).vectorOffset(2) + ((1:n) - 1);
                vectorIndices.vectorIx_RI   = vectorIndices.vectorIx_LI + apertureInfo.totalNumOfLeafPairs;
                vectorIndices.vectorIx_RF   = vectorIndices.vectorIx_LF + apertureInfo.totalNumOfLeafPairs;
            else
                vectorIndices.vectorIx_LI   = updatedInfo.beam(i).shape(j).vectorOffset + ((1:n) - 1);
                vectorIndices.vectorIx_LF   = updatedInfo.beam(i).shape(j).vectorOffset + ((1:n) - 1);
                vectorIndices.vectorIx_RI   = vectorIndices.vectorIx_LI + apertureInfo.totalNumOfLeafPairs;
                vectorIndices.vectorIx_RF   = vectorIndices.vectorIx_LF + apertureInfo.totalNumOfLeafPairs;
            end
        else

            variables.weight_last = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape(j).weight;
            variables.weight_next = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape(j).weight;

            variables.jacobiScale_last    = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape(1).jacobiScale;
            variables.jacobiScale_next    = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape(1).jacobiScale;

            variables.time_last = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).time;
            variables.time_next = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).time;
            variables.time      = updatedInfo.beam(i).time;

            variables.fracFromLastOptI  = fracFromLastOptI;
            variables.fracFromLastOptF  = fracFromLastOptF;
            variables.fracFromNextOptI  = fracFromNextOptI;
            variables.fracFromNextOptF  = fracFromNextOptF;
            variables.fracFromLastOpt   = fracFromLastOpt;

            variables.doseAngleBordersDiff      = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff;
            variables.doseAngleBordersDiff_last = updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).doseAngleBordersDiff;
            variables.doseAngleBordersDiff_next = updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).doseAngleBordersDiff;
            variables.timeFacCurr_last          = updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).timeFacCurr;
            variables.timeFacCurr_next          = updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).timeFacCurr;
            variables.fracFromLastDAO           = updatedInfo.propVMAT.beam(i).fracFromLastDAO;
            variables.timeFracFromLastDAO       = updatedInfo.propVMAT.beam(i).timeFracFromLastDAO;
            variables.timeFracFromNextDAO       = updatedInfo.propVMAT.beam(i).timeFracFromNextDAO;

            vectorIndices.DAOindex_last = updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).DAOIndex;
            vectorIndices.DAOindex_next = updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex;
            vectorIndices.tIx_last      = tIxOffset + vectorIndices.DAOindex_last;
            vectorIndices.tIx_next      = tIxOffset + vectorIndices.DAOindex_next;

            if updatedInfo.continuousAperture
                vectorIndices.vectorIx_LF_last  = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape(j).vectorOffset(2) + ((1:n) - 1);
                vectorIndices.vectorIx_LI_next  = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape(j).vectorOffset(1) + ((1:n) - 1);
                vectorIndices.vectorIx_RF_last  = vectorIndices.vectorIx_LF_last + apertureInfo.totalNumOfLeafPairs;
                vectorIndices.vectorIx_RI_next  = vectorIndices.vectorIx_LI_next + apertureInfo.totalNumOfLeafPairs;
            else
                vectorIndices.vectorIx_LF_last  = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape(j).vectorOffset + ((1:n) - 1);
                vectorIndices.vectorIx_LI_next  = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape(j).vectorOffset + ((1:n) - 1);
                vectorIndices.vectorIx_RF_last  = vectorIndices.vectorIx_LF_last + apertureInfo.totalNumOfLeafPairs;
                vectorIndices.vectorIx_RI_next  = vectorIndices.vectorIx_LI_next + apertureInfo.totalNumOfLeafPairs;
            end
        end

        counters.bixelJApVec_offset = bixelJApVec_offset;

        % calculate bixel weight and derivative in function
        [w, bixelJApVec_vec, bixelJApVec_i, bixelJApVec_j, sumGradSq, shapeMap, counters] = ...
            matRad_bixWeightAndGrad(calcOptions, mlcOptions, variables, vectorIndices, counters, ...
                                    w, bixelJApVec_vec, bixelJApVec_i, bixelJApVec_j, sumGradSq, shapeMap);

        bixelJApVec_offset = counters.bixelJApVec_offset;

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
deleteInd_i = isnan(bixelJApVec_i);
deleteInd_j = bixelJApVec_j == 0;
if ~all(deleteInd_i == deleteInd_j)
    error('Jacobian deletion mismatch');
else
    bixelJApVec_i(deleteInd_i) = [];
    bixelJApVec_j(deleteInd_i) = [];
    bixelJApVec_vec(deleteInd_i) = [];
end
updatedInfo.bixelJApVec = sparse(bixelJApVec_i, bixelJApVec_j, bixelJApVec_vec, numel(apertureInfoVect), updatedInfo.totalNumOfBixels);

end
