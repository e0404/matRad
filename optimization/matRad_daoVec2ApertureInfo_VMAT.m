function updatedInfo = matRad_daoVec2ApertureInfo_VMAT(apertureInfo,apertureInfoVect)
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

w = zeros(apertureInfo.totalNumOfBixels,1);

% initializing variables
updatedInfo = apertureInfo;

updatedInfo.apertureVector = apertureInfoVect;
%if updatedInfo.VMAT
%IandFapertureVector = apertureInfo.IandFapertureVector;
%end


shapeInd = 1;
j = 1; %only 1 shape per beam for VMAT

%indVect = NaN*ones(size(apertureInfoVect));

%change this to eliminate the first unused entries (which pertain to the
%weights of the aprtures, and to make the bixelIndices work when doing VMAT
%(and we need to potentially interpolate between control points)
indVect = NaN*ones(2*apertureInfo.doseTotalNumOfLeafPairs,1);
offset = 0;

% helper function to cope with numerical instabilities through rounding
round2 = @(a,b) round(a*10^b)/10^b;


if ~all([updatedInfo.beam.optimizeBeam])
    for i = 1:numel(updatedInfo.beam)
        if updatedInfo.beam(i).optimizeBeam
            % update the shape weight
            % rescale the weight from the vector using the previous
            % iteration scaling factor
            updatedInfo.beam(i).shape(j).weight = apertureInfoVect(shapeInd)./updatedInfo.beam(i).shape(j).jacobiScale;
            
            updatedInfo.beam(i).MU = updatedInfo.beam(i).shape(j).weight*updatedInfo.weightToMU;
            updatedInfo.beam(i).time = apertureInfoVect(updatedInfo.totalNumOfShapes+updatedInfo.totalNumOfLeafPairs*2+shapeInd)*updatedInfo.beam(i).doseAngleBordersDiff./updatedInfo.beam(i).optAngleBordersDiff;
            updatedInfo.beam(i).gantryRot = updatedInfo.beam(i).doseAngleBordersDiff/updatedInfo.beam(i).time;
            updatedInfo.beam(i).MURate = updatedInfo.beam(i).MU./updatedInfo.beam(i).time;
            
            shapeInd = shapeInd+1;
        end
    end
end

%% update the shapeMaps
% here the new colimator positions are used to create new shapeMaps that
% now include decimal values instead of binary

shapeInd = 1;
% loop over all beams
for i = 1:numel(updatedInfo.beam)
    
    %posOfRightCornerPixel = apertureInfo.beam(i).posOfCornerBixel(1) + (size(apertureInfo.beam(i).bixelIndMap,2)-1)*apertureInfo.bixelWidth;
    
    % pre compute left and right bixel edges
    edges_l = updatedInfo.beam(i).posOfCornerBixel(1)...
        + ([1:size(apertureInfo.beam(i).bixelIndMap,2)]-1-1/2)*updatedInfo.bixelWidth;
    edges_r = updatedInfo.beam(i).posOfCornerBixel(1)...
        + ([1:size(apertureInfo.beam(i).bixelIndMap,2)]-1+1/2)*updatedInfo.bixelWidth;
    
    n = apertureInfo.beam(i).numOfActiveLeafPairs;
    
    % extract leaf position and aperture weight information from vector
    % for optimized beams, this is straightforward
    % for non-optimized beams, interpolate MURate and leaf speed
    if updatedInfo.beam(i).optimizeBeam
        % update the shape weight
        % rescale the weight from the vector using the previous
        % iteration scaling factor
        updatedInfo.beam(i).shape(j).weight = apertureInfoVect(shapeInd)./updatedInfo.beam(i).shape(j).jacobiScale;
        
        updatedInfo.beam(i).MU = updatedInfo.beam(i).shape(j).weight*updatedInfo.weightToMU;
        updatedInfo.beam(i).time = apertureInfoVect(updatedInfo.totalNumOfShapes+updatedInfo.totalNumOfLeafPairs*2+shapeInd)*updatedInfo.beam(i).doseAngleBordersDiff./updatedInfo.beam(i).optAngleBordersDiff;
        updatedInfo.beam(i).gantryRot = updatedInfo.beam(i).doseAngleBordersDiff/updatedInfo.beam(i).time;
        updatedInfo.beam(i).MURate = updatedInfo.beam(i).MU./updatedInfo.beam(i).time;
        
        
        vectorIx_L = updatedInfo.beam(i).shape(j).vectorOffset + ([1:n]-1);
        vectorIx_R = vectorIx_L+apertureInfo.totalNumOfLeafPairs;
        leftLeafPos = apertureInfoVect(vectorIx_L);
        rightLeafPos = apertureInfoVect(vectorIx_R);
        
        
        % update information in shape structure
        updatedInfo.beam(i).shape(j).leftLeafPos  = leftLeafPos;
        updatedInfo.beam(i).shape(j).rightLeafPos = rightLeafPos;
    else
        %MURate is interpolated between MURates of optimized apertures
        updatedInfo.beam(i).MURate = updatedInfo.beam(i).fracFromLastOpt*updatedInfo.beam(updatedInfo.beam(i).lastOptIndex).MURate+(1-updatedInfo.beam(i).fracFromLastOpt)*updatedInfo.beam(updatedInfo.beam(i).nextOptIndex).MURate;
        updatedInfo.beam(i).gantryRot = 1./(updatedInfo.beam(i).timeFracFromLastOpt./updatedInfo.beam(updatedInfo.beam(i).lastOptIndex).gantryRot+updatedInfo.beam(i).timeFracFromNextOpt./updatedInfo.beam(updatedInfo.beam(i).nextOptIndex).gantryRot);
        updatedInfo.beam(i).time = updatedInfo.beam(i).doseAngleBordersDiff./updatedInfo.beam(i).gantryRot;
        
        updatedInfo.beam(i).MU = updatedInfo.beam(i).MURate.*updatedInfo.beam(i).time;
        updatedInfo.beam(i).shape(1).weight = updatedInfo.beam(i).MU./updatedInfo.weightToMU;
        
        vectorIx_L_last = updatedInfo.beam(updatedInfo.beam(i).lastOptIndex).shape(j).vectorOffset + ([1:n]-1);
        vectorIx_R_last = vectorIx_L_last+apertureInfo.totalNumOfLeafPairs;
        leftLeafPos_last = apertureInfoVect(vectorIx_L_last);
        rightLeafPos_last = apertureInfoVect(vectorIx_R_last);
        
        vectorIx_L_next = updatedInfo.beam(updatedInfo.beam(i).nextOptIndex).shape(j).vectorOffset + ([1:n]-1);
        vectorIx_R_next = vectorIx_L_next+apertureInfo.totalNumOfLeafPairs;
        leftLeafPos_next = apertureInfoVect(vectorIx_L_next);
        rightLeafPos_next = apertureInfoVect(vectorIx_R_next);
        
        
        updatedInfo.beam(i).shape(j).leftLeafPos = updatedInfo.beam(i).fracFromLastOpt*leftLeafPos_last+(1-updatedInfo.beam(i).fracFromLastOpt)*leftLeafPos_next;
        updatedInfo.beam(i).shape(j).rightLeafPos = updatedInfo.beam(i).fracFromLastOpt*rightLeafPos_last+(1-updatedInfo.beam(i).fracFromLastOpt)*rightLeafPos_next;
    end
    
    leftLeafPos = updatedInfo.beam(i).shape(j).leftLeafPos;
    rightLeafPos = updatedInfo.beam(i).shape(j).rightLeafPos;
    
    % rounding for numerical stability
    leftLeafPos  = round2(leftLeafPos,10);
    rightLeafPos = round2(rightLeafPos,10);
    
    %%%%POSSIBLY NECESSARY IN CASE LEAF POSITIONS EVER
    %%%%OVERSHOOT MAXIMUMS
    leftLeafPos(leftLeafPos <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(leftLeafPos <= apertureInfo.beam(i).lim_l);
    rightLeafPos(rightLeafPos <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(rightLeafPos <= apertureInfo.beam(i).lim_l);
    leftLeafPos(leftLeafPos >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(leftLeafPos >= apertureInfo.beam(i).lim_r);
    rightLeafPos(rightLeafPos >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(rightLeafPos >= apertureInfo.beam(i).lim_r);
    
    %
    xPosIndLeftLeaf  = round((leftLeafPos - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
    xPosIndRightLeaf = round((rightLeafPos - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
    
    %
    xPosIndLeftLeaf_lim  = floor((apertureInfo.beam(i).lim_l - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth+1);
    xPosIndRightLeaf_lim = ceil((apertureInfo.beam(i).lim_r - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
    
    xPosIndLeftLeaf(xPosIndLeftLeaf <= xPosIndLeftLeaf_lim) = xPosIndLeftLeaf_lim(xPosIndLeftLeaf <= xPosIndLeftLeaf_lim)+1;
    xPosIndRightLeaf(xPosIndRightLeaf >= xPosIndRightLeaf_lim) = xPosIndRightLeaf_lim(xPosIndRightLeaf >= xPosIndRightLeaf_lim)-1;
    
    
    % check limits because of rounding off issues at maximum, i.e.,
    % enforce round(X.5) -> X
    % LeafPos can occasionally go slightly beyond lim_r, so changed
    % == check to >=
    xPosIndLeftLeaf(leftLeafPos >= apertureInfo.beam(i).lim_r) = round(...
        .5 + (leftLeafPos(leftLeafPos >= apertureInfo.beam(i).lim_r) ...
        - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth);
    
    xPosIndRightLeaf(rightLeafPos >= apertureInfo.beam(i).lim_r) = round(...
        .5 + (rightLeafPos(rightLeafPos >= apertureInfo.beam(i).lim_r) ...
        - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth);
    
    %{
            xPosIndLeftLeaf(leftLeafPos == apertureInfo.beam(i).lim_r) = ...
                .5 + (leftLeafPos(leftLeafPos == apertureInfo.beam(i).lim_r) ...
                - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth;
            xPosIndRightLeaf(rightLeafPos == apertureInfo.beam(i).lim_r) = ...
                .5 + (rightLeafPos(rightLeafPos == apertureInfo.beam(i).lim_r) ...
                - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth;
    %}
    
    % find the bixel index that the leaves currently touch
    bixelIndLeftLeaf  = apertureInfo.beam(i).bixelIndMap((xPosIndLeftLeaf-1)*n+[1:n]');
    bixelIndRightLeaf = apertureInfo.beam(i).bixelIndMap((xPosIndRightLeaf-1)*n+[1:n]');
    
    if any(isnan(bixelIndLeftLeaf)) || any(isnan(bixelIndRightLeaf))
        error('cannot map leaf position to bixel index');
    end
    
    % store information in index vector for gradient calculation
    indVect(offset+[1:n]) = bixelIndLeftLeaf;
    indVect(offset+[1:n]+apertureInfo.doseTotalNumOfLeafPairs) = bixelIndRightLeaf;
    offset = offset+n;
    
    % calculate opening fraction for every bixel in shape to construct
    % bixel weight vector
    
    coveredByLeftLeaf  = bsxfun(@minus,leftLeafPos,edges_l)  / updatedInfo.bixelWidth;
    coveredByRightLeaf = bsxfun(@minus,edges_r,rightLeafPos) / updatedInfo.bixelWidth;
    
    tempMap = 1 - (coveredByLeftLeaf  + abs(coveredByLeftLeaf))  / 2 ...
        - (coveredByRightLeaf + abs(coveredByRightLeaf)) / 2;
    
    
    % find open bixels
    tempMapIx = tempMap > 0;
    
    currBixelIx = apertureInfo.beam(i).bixelIndMap(tempMapIx);
    w(currBixelIx) = w(currBixelIx) + tempMap(tempMapIx)*updatedInfo.beam(i).shape(j).weight;
    
    % save the tempMap (we need to apply a positivity operator !)
    updatedInfo.beam(i).shape(j).shapeMap = (tempMap  + abs(tempMap))  / 2;
    
    if updatedInfo.beam(i).optimizeBeam
        
        if apertureInfo.updatePreconditioner
            dijScaleFactor = mean(apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes)./apertureInfo.jacobiScale)/(apertureInfo.bixelWidth);
            
            if apertureInfo.preconditioner
                % "incorrect"
                %updatedInfo.beam(i).shape(j).jacobiScale = sqrt(sum(updatedInfo.beam(i).shape(j).shapeMap(:)));
                % "correct"
                updatedInfo.beam(i).shape(j).jacobiScale = (dijScaleFactor.*apertureInfo.bixelWidth./apertureInfo.beam(i).shape(j).weight).*sqrt(sum(updatedInfo.beam(i).shape(j).shapeMap(:).^2));
            end
            % rescale the vector from the weight using the current
            % iteration scaling factor
            apertureInfoVect(shapeInd) = updatedInfo.beam(i).shape(j).jacobiScale*updatedInfo.beam(i).shape(j).weight;
            
            updatedInfo.jacobiScale(shapeInd) = updatedInfo.beam(i).shape(j).jacobiScale;
        end
        
        % increment shape index
        shapeInd = shapeInd +1;
    end
end


updatedInfo.bixelWeights = w;
updatedInfo.bixelIndices = indVect;
updatedInfo.apertureVector = apertureInfoVect;

%if updatedInfo.VMAT
%updatedInfo.IandFapertureVector = IandFapertureVector;
%end

end