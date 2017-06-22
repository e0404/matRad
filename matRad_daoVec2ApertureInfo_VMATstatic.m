function updatedInfo = matRad_daoVec2ApertureInfo_VMATstatic(apertureInfo,apertureInfoVect,touchingFlag)
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

if nargin < 3
    touchingFlag = 0; %default is 0, it should really only be 1 when called the first time in the leaf sequencing function, never during DAO
end

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

%Only collect this data once, to save time
dimZ = updatedInfo.beam(1).numOfActiveLeafPairs;
numBeams = numel(unique([updatedInfo.beam.gantryAngle]));
leftLeafPoss = nan(dimZ,numBeams); %Each non-interpolated beam should have 1 and only 1 shape
rightLeafPoss = nan(dimZ,numBeams);
optGantryAngles = zeros(1,numBeams);

initBorderGantryAngles = unique([updatedInfo.beam.initAngleBorders]);
initBorderLeftLeafPoss = nan(dimZ,numel(initBorderGantryAngles));

l = 1;
m = 1;
for k = 1:numel(updatedInfo.beam)
    if k ~= 1 && updatedInfo.beam(k).gantryAngle == updatedInfo.beam(k-1).gantryAngle
        continue
    end
    
    if updatedInfo.beam(k).numOfShapes
        leftLeafPoss(:,l) = updatedInfo.beam(k).shape(1).leftLeafPos;
        rightLeafPoss(:,l) = updatedInfo.beam(k).shape(1).rightLeafPos;
        %leftLeafPoss(:,l) = apertureInfoVect(vectorIx);
        %rightLeafPoss(:,l) = apertureInfoVect(vectorIx+updatedInfo.totalNumOfLeafPairs);
        optGantryAngles(l) = updatedInfo.beam(k).gantryAngle;
        
        l = l+1;
    end
    
    %Only important when cleaning up instances of opposing
    %leaves touching.
    if updatedInfo.beam(k).initializeBeam
        if updatedInfo.beam(k).leafDir == 1
            %This means that the current arc sector is moving
            %in the normal direction (L-R).
            initBorderLeftLeafPoss(:,m) = updatedInfo.beam(k).lim_l;
            
        elseif updatedInfo.beam(k).leafDir == -1
            %This means that the current arc sector is moving
            %in the reverse direction (R-L).
            initBorderLeftLeafPoss(:,m) = updatedInfo.beam(k).lim_r;
        end
        m = m+1;
        
        %end of last sector
        if m == numel(initBorderGantryAngles)
            %This gives ending angle of the current sector.
            if updatedInfo.beam(k).leafDir == 1
                %This means that the current arc sector is moving
                %in the normal direction (L-R), so the next arc
                %sector is moving opposite
                initBorderLeftLeafPoss(:,m) = updatedInfo.beam(k).lim_r;
            elseif updatedInfo.beam(k).leafDir == -1
                %This means that the current arc sector is moving
                %in the reverse direction (R-L), so the next
                %arc sector is moving opposite
                initBorderLeftLeafPoss(:,m) = updatedInfo.beam(k).lim_l;
            end
        end
    end
end

%Any time leaf pairs are touching, they are set to
%be in the middle of the field.  Instead, move them
%so that they are still touching, but that they
%follow the motion of the MLCs across the field.
for row = 1:dimZ
    
    touchingInd = find(leftLeafPoss(row,:) == rightLeafPoss(row,:));
    
    if ~exist('leftLeafPossAug','var')
        %leftLeafPossAug = [reshape(mean([leftLeafPoss(:) rightLeafPoss(:)],2),size(leftLeafPoss)),borderLeftLeafPoss];
        leftLeafPossAugTemp = reshape(mean([leftLeafPoss(:) rightLeafPoss(:)],2),size(leftLeafPoss));
        
        numRep = 0;
        repInd = nan(size(optGantryAngles));
        for j = 1:numel(optGantryAngles)
            if any(optGantryAngles(j) == initBorderGantryAngles)
                %replace leaf positions with the ones at
                %the borders (eliminates repetitions)
                numRep = numRep+1;
                %these are the gantry angles that are
                %repeated
                repInd(numRep) = j;
                
                delInd = find(optGantryAngles(j) == initBorderGantryAngles);
                leftLeafPossAugTemp(:,j) = initBorderLeftLeafPoss(:,delInd);
                initBorderLeftLeafPoss(:,delInd) = [];
                initBorderGantryAngles(delInd) = [];
            end
        end
        repInd(isnan(repInd)) = [];
        leftLeafPossAug = [leftLeafPossAugTemp,initBorderLeftLeafPoss];
        gantryAnglesAug = [optGantryAngles,initBorderGantryAngles];
    end
    notTouchingInd = [setdiff(1:numBeams,touchingInd),repInd];
    notTouchingInd = unique(notTouchingInd);
    %make sure to include the repeated ones in the
    %interpolation!
    
    notTouchingIndAug = [notTouchingInd,(1+numel(optGantryAngles)):(numel(optGantryAngles)+numel(initBorderGantryAngles))];
    
    leftLeafPoss(row,touchingInd) = interp1(gantryAnglesAug(notTouchingIndAug),leftLeafPossAug(row,notTouchingIndAug),optGantryAngles(touchingInd));
    rightLeafPoss(row,touchingInd) = leftLeafPoss(row,touchingInd);
end

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
            
            if touchingFlag
                %Perform interpolation
                leftLeafPos = (interp1(optGantryAngles',leftLeafPoss',updatedInfo.beam(i).gantryAngle))';
                rightLeafPos = (interp1(optGantryAngles',rightLeafPoss',updatedInfo.beam(i).gantryAngle))';
                
                %re-update vector in case anything changed from fixing the leaf
                %touching
                vectorIx = updatedInfo.beam(i).shape(j).vectorOffset + ([1:n]-1);
                apertureInfoVect(vectorIx) = leftLeafPos;
                apertureInfoVect(vectorIx+apertureInfo.totalNumOfLeafPairs) = rightLeafPos;
            else
                % extract left and right leaf positions from shape vector
                vectorIx     = updatedInfo.beam(i).shape(j).vectorOffset + ([1:n]-1);
                leftLeafPos  = apertureInfoVect(vectorIx);
                rightLeafPos = apertureInfoVect(vectorIx+apertureInfo.totalNumOfLeafPairs);
            end
            
            % extract left and right leaf positions from shape vector
            vectorIx     = updatedInfo.beam(i).shape(j).vectorOffset + ([1:n]-1);
            leftLeafPos  = apertureInfoVect(vectorIx);
            rightLeafPos = apertureInfoVect(vectorIx+apertureInfo.totalNumOfLeafPairs);
            
            %interpolate initial and final leaf positions
            %if apertureInfo.beam(i).doseAngleOpt(1)
            %initial
            %IandFvectorIx = updatedInfo.beam(i).shape(j).IandFvectorOffset(1) + ([1:n]-1);
            updatedInfo.beam(i).shape(j).leftLeafPos_I = (interp1(optGantryAngles',leftLeafPoss',updatedInfo.beam(i).doseAngleBorders(1)))';
            updatedInfo.beam(i).shape(j).rightLeafPos_I = (interp1(optGantryAngles',rightLeafPoss',updatedInfo.beam(i).doseAngleBorders(1)))';
            
            %IandFapertureVector(IandFvectorIx) = updatedInfo.beam(i).shape(j).leftLeafPos_I;
            %IandFapertureVector(IandFvectorIx+apertureInfo.IandFtotalNumOfLeafPairs) = updatedInfo.beam(i).shape(j).rightLeafPos_I;
            %end
            
            %if apertureInfo.beam(i).doseAngleOpt(2)
            %final
            %IandFvectorIx = updatedInfo.beam(i).shape(j).IandFvectorOffset(2) + ([1:n]-1);
            updatedInfo.beam(i).shape(j).leftLeafPos_F = (interp1(optGantryAngles',leftLeafPoss',updatedInfo.beam(i).doseAngleBorders(2)))';
            updatedInfo.beam(i).shape(j).rightLeafPos_F = (interp1(optGantryAngles',rightLeafPoss',updatedInfo.beam(i).doseAngleBorders(2)))';
            
            %IandFapertureVector(IandFvectorIx) = updatedInfo.beam(i).shape(j).leftLeafPos_F;
            %IandFapertureVector(IandFvectorIx+apertureInfo.IandFtotalNumOfLeafPairs) = updatedInfo.beam(i).shape(j).rightLeafPos_F;
            %end
            
            % update information in shape structure
            updatedInfo.beam(i).shape(j).leftLeafPos  = leftLeafPos;
            updatedInfo.beam(i).shape(j).rightLeafPos = rightLeafPos;
            
            
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