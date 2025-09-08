function updatedInfo = matRad_daoVec2ApertureInfo_VMATrecalcDynamic(apertureInfo,apertureInfoVect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to translate the vector representation of the aperture
% shape and weight into an aperture info struct. At the same time, the
% updated bixel weight vector w is computed and a vector listing the
% correspondence between leaf tips and bixel indices for gradient
% calculation
%
% call
%   updatedInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVect)
%
% input
%   apertureInfo:     aperture shape info struct
%   apertureInfoVect: aperture weights and shapes parameterized as vector
%   touchingFlag:     if this is one, clean up instances of leaf touching,
%                     otherwise, do not
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

% function to update the apertureInfo struct after the each iteraton of the
% optimization

% initializing variables
updatedInfo = apertureInfo;

updatedInfo.apertureVector = apertureInfoVect;

% options for bixel and Jacobian calculation
mlcOptions.bixelWidth = apertureInfo.bixelWidth;
calcOptions.continuousAperture = updatedInfo.propVMAT.continuousAperture;
vectorIndices.totalNumOfShapes = apertureInfo.totalNumOfShapes;

w = cell(apertureInfo.numPhases,1);
w(:) = {zeros(apertureInfo.totalNumOfBixels,1)};

if updatedInfo.runVMAT && ~all([updatedInfo.propVMAT.beam.DAOBeam])
    j = 1;
    for i = 1:numel(updatedInfo.beam)
        if updatedInfo.propVMAT.beam(i).DAOBeam
            % update the shape weight
            % rescale the weight from the vector using the previous
            % iteration scaling factor
            updatedInfo.beam(i).shape(j).weight = apertureInfoVect(updatedInfo.beam(i).shape(j).weightOffset)./updatedInfo.beam(i).shape(j).jacobiScale;
            
            updatedInfo.beam(i).shape(j).MU = updatedInfo.beam(i).shape(j).weight*updatedInfo.weightToMU;
            updatedInfo.beam(i).time = apertureInfoVect((updatedInfo.totalNumOfShapes+updatedInfo.totalNumOfLeafPairs*2)+updatedInfo.propVMAT.beam(i).DAOIndex)*updatedInfo.propVMAT.beam(i).timeFacCurr;
            updatedInfo.beam(i).gantryRot = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff/updatedInfo.beam(i).time;
            updatedInfo.beam(i).shape(j).MURate = updatedInfo.beam(i).shape(j).MU./updatedInfo.beam(i).time;
        end
    end
end

bixelJApVec_vec = cell(apertureInfo.numPhases,1);

% dummy variables
bixelJApVec_i = 0;
bixelJApVec_j = 0;
bixelJApVec_offset = 0;
counters.bixelJApVec_offset = bixelJApVec_offset;

% Interpolate segment between adjacent optimized gantry angles.
% Include in updatedInfo, but NOT the vector (since these are not
% optimized by DAO).  Also update bixel weights to include these.


%% update the shapeMaps
% here the new colimator positions are used to create new shapeMaps that
% now include decimal values instead of binary

calcOptions.saveJacobian = false;

% loop over all beams
for i = 1:numel(updatedInfo.beam)
    
    %posOfRightCornerPixel = apertureInfo.beam(i).posOfCornerBixel(1) + (size(apertureInfo.beam(i).bixelIndMap,2)-1)*apertureInfo.bixelWidth;
    
    % pre compute left and right bixel edges
    edges_l = updatedInfo.beam(i).posOfCornerBixel(1)...
        + ([1:size(apertureInfo.beam(i).bixelIndMap,2)]-1-1/2)*updatedInfo.bixelWidth;
    edges_r = updatedInfo.beam(i).posOfCornerBixel(1)...
        + ([1:size(apertureInfo.beam(i).bixelIndMap,2)]-1+1/2)*updatedInfo.bixelWidth;
    
    % get dimensions of 2d matrices that store shape/bixel information
    n = apertureInfo.beam(i).numOfActiveLeafPairs;
    
    %weightFactor_I = updatedInfo.propVMAT.beam(i).doseAngleBorderCentreDiff(1)./updatedInfo.propVMAT.beam(i).doseAngleBordersDiff;
    %weightFactor_F = updatedInfo.propVMAT.beam(i).doseAngleBorderCentreDiff(2)./updatedInfo.propVMAT.beam(i).doseAngleBordersDiff;
    
    % we are necessarily doing VMAT
    numOfShapes = 1;
    
    mlcOptions.lim_l = apertureInfo.beam(i).lim_l;
    mlcOptions.lim_r = apertureInfo.beam(i).lim_r;
    mlcOptions.edges_l = edges_l;
    mlcOptions.edges_r = edges_r;
    mlcOptions.centres = (edges_l+edges_r)/2;
    mlcOptions.widths = edges_r-edges_l;
    mlcOptions.n = n;
    mlcOptions.numBix = size(apertureInfo.beam(i).bixelIndMap,2);
    mlcOptions.bixelIndMap = apertureInfo.beam(i).bixelIndMap;
    calcOptions.DAOBeam = updatedInfo.propVMAT.beam(i).DAOBeam;
    
    % loop over all shapes
    for j = 1:numOfShapes
        
        % shapeMap
        shapeMap_I = zeros(size(updatedInfo.beam(i).bixelIndMap));
        shapeMap_F = zeros(size(updatedInfo.beam(i).bixelIndMap));
        % sumGradSq
        sumGradSq = 0;
        
        % no need to update weights or anything from the vector, just
        % extract the weights and leaf positions from the apertureInfo
        
        weight = updatedInfo.beam(i).shape(j).weight;
        if isfield(updatedInfo.beam(i).shape(j),'weight_I')
            weight_I = updatedInfo.beam(i).shape(j).weight_I;
            weight_F = updatedInfo.beam(i).shape(j).weight_F;
        else
            %only happens at original angular resolution
            weight_I = weight.*updatedInfo.beam(i).doseAngleBorderCentreDiff(1)./updatedInfo.beam(i).doseAngleBordersDiff;
            weight_F = weight.*updatedInfo.beam(i).doseAngleBorderCentreDiff(2)./updatedInfo.beam(i).doseAngleBordersDiff;
        end
        
        if weight_I+weight_F ~= weight
            %sometimes the sum is different than one by ~10^-16
            %(rounding error in the division)
            weight_F = weight-weight_I;
        end
        
        %% enter in variables and options
        
        %%%%%%%%%%%%%%%%
        %do initial and final arc separately, more accurate
        %calculation
        
        %INITIAL
        variables.weight_I          = weight_I;
        variables.weight_F          = weight_I;
        variables.weightFactor_I    = 1/2;
        variables.weightFactor_F    = 1/2;
        
        if updatedInfo.propVMAT.continuousAperture
            variables.leftLeafPos_I     = updatedInfo.beam(i).shape(j).leftLeafPos_I;
            variables.leftLeafPos_F     = updatedInfo.beam(i).shape(j).leftLeafPos;
            variables.rightLeafPos_I    = updatedInfo.beam(i).shape(j).rightLeafPos_I;
            variables.rightLeafPos_F    = updatedInfo.beam(i).shape(j).rightLeafPos;
        else
            variables.leftLeafPos_I     = updatedInfo.beam(i).shape(j).leftLeafPos;
            variables.leftLeafPos_F     = updatedInfo.beam(i).shape(j).leftLeafPos;
            variables.rightLeafPos_I    = updatedInfo.beam(i).shape(j).rightLeafPos;
            variables.rightLeafPos_F    = updatedInfo.beam(i).shape(j).rightLeafPos;
        end
        
        % calculate bixel weight and derivative in function
        [w,~,bixelJApVec_i,bixelJApVec_j,sumGradSq,shapeMap_I,counters] = ...
            matRad_bixWeightAndGrad(calcOptions,mlcOptions,variables,vectorIndices,counters,w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j,sumGradSq,shapeMap_I);
        
        %FINAL
        variables.weight_I          = weight_F;
        variables.weight_F          = weight_F;
        variables.weightFactor_I    = 1/2;
        variables.weightFactor_F    = 1/2;
        
        if updatedInfo.propVMAT.continuousAperture
            variables.leftLeafPos_I     = updatedInfo.beam(i).shape(j).leftLeafPos;
            variables.leftLeafPos_F     = updatedInfo.beam(i).shape(j).leftLeafPos_F;
            variables.rightLeafPos_I    = updatedInfo.beam(i).shape(j).rightLeafPos;
            variables.rightLeafPos_F    = updatedInfo.beam(i).shape(j).rightLeafPos_F;
        else
            variables.leftLeafPos_I     = updatedInfo.beam(i).shape(j).leftLeafPos;
            variables.leftLeafPos_F     = updatedInfo.beam(i).shape(j).leftLeafPos;
            variables.rightLeafPos_I    = updatedInfo.beam(i).shape(j).rightLeafPos;
            variables.rightLeafPos_F    = updatedInfo.beam(i).shape(j).rightLeafPos;
        end
        
        % calculate bixel weight and derivative in function
        [w,~,bixelJApVec_i,bixelJApVec_j,sumGradSq,shapeMap_F,counters] = ...
            matRad_bixWeightAndGrad(calcOptions,mlcOptions,variables,vectorIndices,counters,w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j,sumGradSq,shapeMap_F);
        
        % save the tempMap
        shapeMap = shapeMap_I+shapeMap_F;
        updatedInfo.beam(i).shape(j).shapeMap = shapeMap;
    end
    
end

% save bixelWeight, apertureVector, and Jacobian between the two
updatedInfo.bixelWeights = w;
updatedInfo.apertureVector = apertureInfoVect;

end