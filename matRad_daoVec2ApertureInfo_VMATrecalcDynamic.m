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

w = zeros(apertureInfo.totalNumOfBixels,1);

% initializing variables
updatedInfo = apertureInfo;

updatedInfo.apertureVector = apertureInfoVect;
%if updatedInfo.VMAT
%IandFapertureVector = apertureInfo.IandFapertureVector;
%end


shapeInd = 1;

%indVect = NaN*ones(size(apertureInfoVect));

%change this to eliminate the first unused entries (which pertain to the
%weights of the aprtures, and to make the bixelIndices work when doing VMAT
%(and we need to potentially interpolate between control points)
indVect = NaN*ones(2*apertureInfo.doseTotalNumOfLeafPairs,1);
offset = 0;

% helper function to cope with numerical instabilities through rounding
round2 = @(a,b) round(a*10^b)/10^b;


%Interpolate segment between adjacent optimized gantry angles.
% Include in updatedInfo, but NOT the vector (since these are not
% optimized by DAO).  Also update bixel weights to include these.


%% update the shapeMaps
% here the new colimator positions are used to create new shapeMaps that
% now include decimal values instead of binary

% loop over all beams
for i = 1:numel(updatedInfo.beam)
    
    %posOfRightCornerPixel = apertureInfo.beam(i).posOfCornerBixel(1) + (size(apertureInfo.beam(i).bixelIndMap,2)-1)*apertureInfo.bixelWidth;
    
    % pre compute left and right bixel edges
    edges_l = updatedInfo.beam(i).posOfCornerBixel(1)...
        + ([1:size(apertureInfo.beam(i).bixelIndMap,2)]-1-1/2)*updatedInfo.bixelWidth;
    edges_r = updatedInfo.beam(i).posOfCornerBixel(1)...
        + ([1:size(apertureInfo.beam(i).bixelIndMap,2)]-1+1/2)*updatedInfo.bixelWidth;
    
    if updatedInfo.beam(i).numOfShapes ~= 0 %numOfShapes is 0 for interpolated beams!
        
        % loop over all shapes
        for j = 1:updatedInfo.beam(i).numOfShapes
            % update the shape weight
            updatedInfo.beam(i).shape(j).weight = apertureInfoVect(shapeInd);
            
            updatedInfo.beam(i).MU = updatedInfo.beam(i).shape(j).weight*updatedInfo.weightToMU;
            updatedInfo.beam(i).time = apertureInfoVect(updatedInfo.totalNumOfShapes+updatedInfo.totalNumOfLeafPairs*2+shapeInd);
            updatedInfo.beam(i).gantryRot = updatedInfo.beam(i).optAngleBordersDiff/updatedInfo.beam(i).time;
            updatedInfo.beam(i).MURate = updatedInfo.beam(i).MU*updatedInfo.beam(i).gantryRot/updatedInfo.beam(i).doseAngleBordersDiff;
            
            % get dimensions of 2d matrices that store shape/bixel information
            n = apertureInfo.beam(i).numOfActiveLeafPairs;
            
            leftLeafPosI_I = min([updatedInfo.beam(i).shape(j).leftLeafPos_I,updatedInfo.beam(i).shape(j).leftLeafPos],[],2);
            leftLeafPosF_I = max([updatedInfo.beam(i).shape(j).leftLeafPos_I,updatedInfo.beam(i).shape(j).leftLeafPos],[],2);
            rightLeafPosI_I = min([updatedInfo.beam(i).shape(j).rightLeafPos_I,updatedInfo.beam(i).shape(j).rightLeafPos],[],2);
            rightLeafPosF_I = max([updatedInfo.beam(i).shape(j).rightLeafPos_I,updatedInfo.beam(i).shape(j).rightLeafPos],[],2);
            leftLeafPosI_F = min([updatedInfo.beam(i).shape(j).leftLeafPos,updatedInfo.beam(i).shape(j).leftLeafPos_F],[],2);
            leftLeafPosF_F = max([updatedInfo.beam(i).shape(j).leftLeafPos,updatedInfo.beam(i).shape(j).leftLeafPos_F],[],2);
            rightLeafPosI_F = min([updatedInfo.beam(i).shape(j).rightLeafPos,updatedInfo.beam(i).shape(j).rightLeafPos_F],[],2);
            rightLeafPosF_F = max([updatedInfo.beam(i).shape(j).rightLeafPos,updatedInfo.beam(i).shape(j).rightLeafPos_F],[],2);
            
            %{
            leftLeafPosI_I = round2(leftLeafPosI_I,10);
            leftLeafPosF_I = round2(leftLeafPosF_I,10);
            rightLeafPosI_I = round2(rightLeafPosI_I,10);
            rightLeafPosF_I = round2(rightLeafPosF_I,10);
            leftLeafPosI_F = round2(leftLeafPosI_F,10);
            leftLeafPosF_F = round2(leftLeafPosF_F,10);
            rightLeafPosI_F = round2(rightLeafPosI_F,10);
            rightLeafPosF_F = round2(rightLeafPosF_F,10);
            %}
            
            leftLeafPosI_I(leftLeafPosI_I <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(leftLeafPosI_I <= apertureInfo.beam(i).lim_l);
            leftLeafPosF_I(leftLeafPosF_I <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(leftLeafPosF_I <= apertureInfo.beam(i).lim_l);
            rightLeafPosI_I(rightLeafPosI_I <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(rightLeafPosI_I <= apertureInfo.beam(i).lim_l);
            rightLeafPosF_I(rightLeafPosF_I <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(rightLeafPosF_I <= apertureInfo.beam(i).lim_l);
            leftLeafPosI_I(leftLeafPosI_I >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(leftLeafPosI_I >= apertureInfo.beam(i).lim_r);
            leftLeafPosF_I(leftLeafPosF_I >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(leftLeafPosF_I >= apertureInfo.beam(i).lim_r);
            rightLeafPosI_I(rightLeafPosI_I >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(rightLeafPosI_I >= apertureInfo.beam(i).lim_r);
            rightLeafPosF_I(rightLeafPosF_I >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(rightLeafPosF_I >= apertureInfo.beam(i).lim_r);
            leftLeafPosI_F(leftLeafPosI_F <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(leftLeafPosI_F <= apertureInfo.beam(i).lim_l);
            leftLeafPosF_F(leftLeafPosF_F <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(leftLeafPosF_F <= apertureInfo.beam(i).lim_l);
            rightLeafPosI_F(rightLeafPosI_F <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(rightLeafPosI_F <= apertureInfo.beam(i).lim_l);
            rightLeafPosF_F(rightLeafPosF_F <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(rightLeafPosF_F <= apertureInfo.beam(i).lim_l);
            leftLeafPosI_F(leftLeafPosI_F >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(leftLeafPosI_F >= apertureInfo.beam(i).lim_r);
            leftLeafPosF_F(leftLeafPosF_F >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(leftLeafPosF_F >= apertureInfo.beam(i).lim_r);
            rightLeafPosI_F(rightLeafPosI_F >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(rightLeafPosI_F >= apertureInfo.beam(i).lim_r);
            rightLeafPosF_F(rightLeafPosF_F >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(rightLeafPosF_F >= apertureInfo.beam(i).lim_r);
            
            
            %store these for gradient?
            xPosIndLeftLeafI_I  = round((leftLeafPosI_I - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
            xPosIndLeftLeafF_I  = round((leftLeafPosF_I - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
            xPosIndRightLeafI_I = round((rightLeafPosI_I - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
            xPosIndRightLeafF_I = round((rightLeafPosF_I - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
            xPosIndLeftLeafI_F  = round((leftLeafPosI_F - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
            xPosIndLeftLeafF_F  = round((leftLeafPosF_F - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
            xPosIndRightLeafI_F = round((rightLeafPosI_F - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
            xPosIndRightLeafF_F = round((rightLeafPosF_F - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
            
            xPosLinearIndLeftLeafI_I = sub2ind([n size(edges_l)],(1:n)',xPosIndLeftLeafI_I);
            xPosLinearIndLeftLeafF_I = sub2ind([n size(edges_l)],(1:n)',xPosIndLeftLeafF_I);
            xPosLinearIndRightLeafI_I = sub2ind([n size(edges_l)],(1:n)',xPosIndRightLeafI_I);
            xPosLinearIndRightLeafF_I = sub2ind([n size(edges_l)],(1:n)',xPosIndRightLeafF_I);
            xPosLinearIndLeftLeafI_F = sub2ind([n size(edges_l)],(1:n)',xPosIndLeftLeafI_F);
            xPosLinearIndLeftLeafF_F = sub2ind([n size(edges_l)],(1:n)',xPosIndLeftLeafF_F);
            xPosLinearIndRightLeafI_F = sub2ind([n size(edges_l)],(1:n)',xPosIndRightLeafI_F);
            xPosLinearIndRightLeafF_F = sub2ind([n size(edges_l)],(1:n)',xPosIndRightLeafF_F);
            
            %%%%%%%%%%%%%%%%
            %do initial and final arc separately, more accurate
            %calculation?
            
            %INITIAL
            %calculate fraction of fluence uncovered by left leaf
            %initial computation
            uncoveredByLeftLeaf_I = bsxfun(@minus,(edges_l+edges_r)/2,leftLeafPosI_I)./repmat(leftLeafPosF_I-leftLeafPosI_I,1,size(edges_l,2));
            %correct for overshoot in initial and final leaf positions
            uncoveredByLeftLeaf_I(xPosLinearIndLeftLeafI_I) = uncoveredByLeftLeaf_I(xPosLinearIndLeftLeafI_I)+(leftLeafPosI_I-edges_l(xPosIndLeftLeafI_I)').^2./((leftLeafPosF_I-leftLeafPosI_I).*(edges_r(xPosIndLeftLeafI_I)'-edges_l(xPosIndLeftLeafI_I)').*2);
            uncoveredByLeftLeaf_I(xPosLinearIndLeftLeafF_I) = uncoveredByLeftLeaf_I(xPosLinearIndLeftLeafF_I)-(edges_r(xPosIndLeftLeafF_I)'-leftLeafPosF_I).^2./((leftLeafPosF_I-leftLeafPosI_I).*(edges_r(xPosIndLeftLeafF_I)'-edges_l(xPosIndLeftLeafF_I)').*2);
            %round <0 to 0, >1 to 1
            uncoveredByLeftLeaf_I(uncoveredByLeftLeaf_I < 0) = 0;
            uncoveredByLeftLeaf_I(uncoveredByLeftLeaf_I > 1) = 1;
            %uncoveredByLeftLeaf_I = (uncoveredByLeftLeaf_I+abs(uncoveredByLeftLeaf_I))/2;
            %uncoveredByLeftLeaf_I = (uncoveredByLeftLeaf_I+1-abs(uncoveredByLeftLeaf_I-1))/2;
            
            %calculate fraction of fluence covered by right leaf
            %initial computation
            coveredByRightLeaf_I = bsxfun(@minus,(edges_l+edges_r)/2,rightLeafPosI_I)./repmat(rightLeafPosF_I-rightLeafPosI_I,1,size(edges_l,2));
            %correct for overshoot in initial and final leaf positions
            coveredByRightLeaf_I(xPosLinearIndRightLeafI_I) = coveredByRightLeaf_I(xPosLinearIndRightLeafI_I)+(rightLeafPosI_I-edges_l(xPosIndRightLeafI_I)').^2./((rightLeafPosF_I-rightLeafPosI_I).*(edges_r(xPosIndRightLeafI_I)'-edges_l(xPosIndRightLeafI_I)').*2);
            coveredByRightLeaf_I(xPosLinearIndRightLeafF_I) = coveredByRightLeaf_I(xPosLinearIndRightLeafF_I)-(edges_r(xPosIndRightLeafF_I)'-rightLeafPosF_I).^2./((rightLeafPosF_I-rightLeafPosI_I).*(edges_r(xPosIndRightLeafF_I)'-edges_l(xPosIndRightLeafF_I)').*2);
            %round <0 to 0, >1 to 1
            coveredByRightLeaf_I(coveredByRightLeaf_I < 0) = 0;
            coveredByRightLeaf_I(coveredByRightLeaf_I > 1) = 1;
            
            %FINAL
            %calculate fraction of fluence uncovered by left leaf
            %initial computation
            uncoveredByLeftLeaf_F = bsxfun(@minus,(edges_l+edges_r)/2,leftLeafPosI_F)./repmat(leftLeafPosF_F-leftLeafPosI_F,1,size(edges_l,2));
            %correct for overshoot in initial and final leaf positions
            uncoveredByLeftLeaf_F(xPosLinearIndLeftLeafI_F) = uncoveredByLeftLeaf_F(xPosLinearIndLeftLeafI_F)+(leftLeafPosI_F-edges_l(xPosIndLeftLeafI_F)').^2./((leftLeafPosF_F-leftLeafPosI_F).*(edges_r(xPosIndLeftLeafI_F)'-edges_l(xPosIndLeftLeafI_F)').*2);
            uncoveredByLeftLeaf_F(xPosLinearIndLeftLeafF_F) = uncoveredByLeftLeaf_F(xPosLinearIndLeftLeafF_F)-(edges_r(xPosIndLeftLeafF_F)'-leftLeafPosF_F).^2./((leftLeafPosF_F-leftLeafPosI_F).*(edges_r(xPosIndLeftLeafF_F)'-edges_l(xPosIndLeftLeafF_F)').*2);
            %round <0 to 0, >1 to 1
            uncoveredByLeftLeaf_F(uncoveredByLeftLeaf_F < 0) = 0;
            uncoveredByLeftLeaf_F(uncoveredByLeftLeaf_F > 1) = 1;
            
            %calculate fraction of fluence covered by right leaf
            %initial computation
            coveredByRightLeaf_F = bsxfun(@minus,(edges_l+edges_r)/2,rightLeafPosI_F)./repmat(rightLeafPosF_F-rightLeafPosI_F,1,size(edges_l,2));
            %correct for overshoot in initial and final leaf positions
            coveredByRightLeaf_F(xPosLinearIndRightLeafI_F) = coveredByRightLeaf_F(xPosLinearIndRightLeafI_F)+(rightLeafPosI_F-edges_l(xPosIndRightLeafI_F)').^2./((rightLeafPosF_F-rightLeafPosI_F).*(edges_r(xPosIndRightLeafI_F)'-edges_l(xPosIndRightLeafI_F)').*2);
            coveredByRightLeaf_F(xPosLinearIndRightLeafF_F) = coveredByRightLeaf_F(xPosLinearIndRightLeafF_F)-(edges_r(xPosIndRightLeafF_F)'-rightLeafPosF_F).^2./((rightLeafPosF_F-rightLeafPosI_F).*(edges_r(xPosIndRightLeafF_F)'-edges_l(xPosIndRightLeafF_F)').*2);
            %round <0 to 0, >1 to 1
            coveredByRightLeaf_F(coveredByRightLeaf_F < 0) = 0;
            coveredByRightLeaf_F(coveredByRightLeaf_F > 1) = 1;
            
            %fluence is equal to fluence not covered by left leaf minus
            %fluence covered by left leaf
            %round for numerical stability (some bixels which should be
            %0 get arbitrarily close to, but not equal to 0)
            
            tempMap_I = uncoveredByLeftLeaf_I-coveredByRightLeaf_I;
            tempMap_F = uncoveredByLeftLeaf_F-coveredByRightLeaf_F;
            tempMap_I(isnan(tempMap_I)) = 0;
            tempMap_F(isnan(tempMap_F)) = 0;
            
            if isfield(updatedInfo.beam(i).shape(j),'weight_I')
                weight_I = updatedInfo.beam(i).shape(j).weight_I./updatedInfo.beam(i).shape(j).weight;
                weight_F = updatedInfo.beam(i).shape(j).weight_F./updatedInfo.beam(i).shape(j).weight;
            else
                %only happens at original angular resolution
                weight_I = updatedInfo.beam(i).doseAngleBorderCentreDiff(1)./updatedInfo.beam(i).doseAngleBordersDiff;
                weight_F = updatedInfo.beam(i).doseAngleBorderCentreDiff(2)./updatedInfo.beam(i).doseAngleBordersDiff;
            end
            
            if weight_I+weight_F ~= 1
                %sometimes the sum is different than one by ~10^-16
                %(rounding error in the division)
                weight_F = 1-weight_I;
            end
            
            tempMap = weight_I.*tempMap_I+weight_F.*tempMap_F;
            
            if any(tempMap(:) < 0 | tempMap(:) > 1)
                error('map exceeds regular thresholds of 0 and 1');
            end
            
            % find open bixels
            tempMapIx = tempMap > 0;
            
            currBixelIx = apertureInfo.beam(i).bixelIndMap(tempMapIx);
            w(currBixelIx) = w(currBixelIx) + tempMap(tempMapIx)*updatedInfo.beam(i).shape(j).weight;
            
            % save the tempMap (we need to apply a positivity operator !)
            updatedInfo.beam(i).shape(j).shapeMap = tempMap;
            %updatedInfo.beam(i).shape(j).shapeMap = (tempMap  + abs(tempMap))  / 2;
            
            % increment shape index
            shapeInd = shapeInd +1;
        end
        
    else
        %This should only occur for VMAT subchildren angles, i.e., not
        %independently optimized
        
        % get dimensions of 2d matrices that store shape/bixel information
        n = apertureInfo.beam(i).numOfActiveLeafPairs;
        
        %Perform interpolation
        leftLeafPos = (interp1(optGantryAngles',leftLeafPoss',updatedInfo.beam(i).gantryAngle))';
        rightLeafPos = (interp1(optGantryAngles',rightLeafPoss',updatedInfo.beam(i).gantryAngle))';
        
        %MURate is interpolated between MURates of optimized apertures
        updatedInfo.beam(i).MURate = updatedInfo.beam(i).fracFromLastOpt*updatedInfo.beam(updatedInfo.beam(i).lastOptIndex).MURate+(1-updatedInfo.beam(i).fracFromLastOpt)*updatedInfo.beam(updatedInfo.beam(i).nextOptIndex).MURate;
        updatedInfo.beam(i).gantryRot = updatedInfo.beam(i).fracFromLastOpt*updatedInfo.beam(updatedInfo.beam(i).lastOptIndex).gantryRot+(1-updatedInfo.beam(i).fracFromLastOpt)*updatedInfo.beam(updatedInfo.beam(i).nextOptIndex).gantryRot;
        
        updatedInfo.beam(i).MU = updatedInfo.beam(i).MURate*updatedInfo.beam(i).doseAngleBordersDiff/updatedInfo.beam(i).gantryRot;
        updatedInfo.beam(i).shape(1).weight = updatedInfo.beam(i).MU./updatedInfo.weightToMU;
        
        % update information in shape structure
        updatedInfo.beam(i).shape(1).leftLeafPos  = leftLeafPos;
        updatedInfo.beam(i).shape(1).rightLeafPos = rightLeafPos;
        
        %The following is taken from the non-VMAT case (j->1, since there is only 1
        %shape per beam in VMAT)
        % rounding for numerical stability
        leftLeafPos  = round2(leftLeafPos,10);
        rightLeafPos = round2(rightLeafPos,10);
        
        %
        xPosIndLeftLeaf  = round((leftLeafPos - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
        xPosIndRightLeaf = round((rightLeafPos - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
        
        % check limits because of rounding off issues at maximum, i.e.,
        % enforce round(X.5) -> X
        xPosIndLeftLeaf(leftLeafPos == apertureInfo.beam(i).lim_r) = ...
            .5 + (leftLeafPos(leftLeafPos == apertureInfo.beam(i).lim_r) ...
            - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth;
        xPosIndRightLeaf(rightLeafPos == apertureInfo.beam(i).lim_r) = ...
            .5 + (rightLeafPos(rightLeafPos == apertureInfo.beam(i).lim_r) ...
            - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth;
        
        % find the bixel index that the leaves currently touch
        bixelIndLeftLeaf  = apertureInfo.beam(i).bixelIndMap((xPosIndLeftLeaf-1)*dimZ+[1:dimZ]');
        bixelIndRightLeaf = apertureInfo.beam(i).bixelIndMap((xPosIndRightLeaf-1)*dimZ+[1:dimZ]');
        
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
        w(currBixelIx) = w(currBixelIx) + tempMap(tempMapIx)*updatedInfo.beam(i).shape(1).weight;
        
        % save the tempMap (we need to apply a positivity operator !)
        updatedInfo.beam(i).shape(1).shapeMap = (tempMap  + abs(tempMap))  / 2;
        
    end
    
    
end

updatedInfo.bixelWeights = w;
updatedInfo.bixelIndices = indVect;
updatedInfo.apertureVector = apertureInfoVect;

%if updatedInfo.VMAT
%updatedInfo.IandFapertureVector = IandFapertureVector;
%end

end