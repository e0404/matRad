function [updatedInfo,w,indVect] = matRad_vec2ShapeInfo(shapeInfo,shapeInfoVect)

% function to update the shapeInfo struct after the each iteraton of the
% optimization

w = zeros(shapeInfo.totalNumOfBixels,1);

% initializing variables
updatedInfo = shapeInfo;

shapeInd = 1;

indVect = NaN*ones(shapeInfo.totalNumOfShapes + shapeInfo.totalNumOfLeafPairs,1);

%% update the shapeMaps
% here the new colimator positions are used to create new shapeMaps that
% now include decimal values instead of binary

% loop over all beams
for i=1:numel(updatedInfo.beam)
    
    posOfRightCornerPixel = shapeInfo.beam(i).posOfCornerBixel(1) + (size(shapeInfo.beam(1).bixelIndMap,2)-1)*shapeInfo.bixelWidth;
    
    % loop over all shapes
    for j=1:updatedInfo.beam(i).numOfShapes
        
        % update the shape weight
        updatedInfo.beam(i).shape(j).weight = shapeInfoVect(shapeInd);
        
        % create a new shapeMap
        tempMap = updatedInfo.beam(i).shape(j).shapeMap;
        
        for k = 1:updatedInfo.beam(i).numOfActiveLeafPairs
            % fill the new shapeMap
            vectorIx = updatedInfo.beam(i).shape(j).vectorOffset + (k-1);
            leftLeafPos = shapeInfoVect(vectorIx);
            rightLeafPos = shapeInfoVect(vectorIx+shapeInfo.totalNumOfLeafPairs);

            updatedInfo.beam(i).shape(j).leftLeafPos(k) = leftLeafPos;
            updatedInfo.beam(i).shape(j).rightLeafPos(k) = rightLeafPos;
            
            xPosIndLeftLeaf  = size(shapeInfo.beam(1).bixelIndMap,2) - round((posOfRightCornerPixel-leftLeafPos)/shapeInfo.bixelWidth);
            xPosIndRightLeaf = round((rightLeafPos - shapeInfo.beam(i).posOfCornerBixel(1))/shapeInfo.bixelWidth) + 1;
            
            % check limits
            if leftLeafPos == shapeInfo.beam(i).lim_l(k) % if leaf at min_let
                % use first possible bixel
                xPosIndLeftLeaf = find(~isnan(shapeInfo.beam(i).bixelIndMap(k,:)),1,'first');
            end
            if rightLeafPos == shapeInfo.beam(i).lim_r(k) % if leaf at min_let
                % use last possible bixel
                xPosIndRightLeaf = find(~isnan(shapeInfo.beam(i).bixelIndMap(k,:)),1,'last');
            end
            
            bixelIndLeftLeaf  = shapeInfo.beam(i).bixelIndMap(k,xPosIndLeftLeaf);
            bixelIndRightLeaf = shapeInfo.beam(i).bixelIndMap(k,xPosIndRightLeaf);
            
            indVect(shapeInfo.beam(i).shape(j).vectorOffset+k-1) = bixelIndLeftLeaf;
            indVect(shapeInfo.beam(i).shape(j).vectorOffset+k-1+shapeInfo.totalNumOfLeafPairs) = bixelIndRightLeaf;
            
            % left bixel edge
            edges_l = updatedInfo.beam(i).posOfCornerBixel(1)...
                          +([1:size(tempMap,2)]-1-1/2)*updatedInfo.bixelWidth;
            % right bixel edge
            edges_r = updatedInfo.beam(i).posOfCornerBixel(1)...
                        +([1:size(tempMap,2)]-1+1/2)*updatedInfo.bixelWidth;
    
            % determine the necessary left/right position for the
            % opening fraction. This will only correspond to the leaf
            % positions if the current bixel is partially covered by
            % the current leaf pair
            x_l = max([edges_l;min(leftLeafPos,edges_r)]);
            x_r = min([edges_r;max(rightLeafPos,edges_l)]);
                
            openingFrac = (x_r-x_l)/updatedInfo.bixelWidth;
            tempMap(k,:) = openingFrac;
            
            currBixelIx = shapeInfo.beam(i).bixelIndMap(k,openingFrac > 0);
            w(currBixelIx) = w(currBixelIx) + openingFrac(openingFrac > 0)'*updatedInfo.beam(i).shape(j).weight;
            
        end
        
        % save the tempMap
        updatedInfo.beam(i).shape(j).shapeMap = tempMap;
        
        % increment shape index
        shapeInd = shapeInd +1;
    end
    
end

end