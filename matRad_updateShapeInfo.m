function [updatedInfo] = matRad_updateShapeInfo(shapeInfo,shapeInfoVect)

% function to update the shapeInfo struct after the each iteraton of the
% optimization

% initializing variables
updatedInfo = shapeInfo;
totalNumOfLeafPairs = [shapeInfo.beam(:).numOfShapes]*...
                            [shapeInfo.beam(:).numOfActiveLeafPairs]';
shapeInd = 1;

%% update the shapeMaps
% here the new colimator positions are used to create new shapeMaps that
% now include decimal values instead of binary

% loop over all beams
for i=1:numel(updatedInfo.beam)
    % loop over all shapes
    for j=1:updatedInfo.beam(i).numOfShapes
        % create a new shapeMap
        tempMap = updatedInfo.beam(i).shape(j).shapeMap;
        % get beam and shape Ix
%         beamIx = addInfoVect(:,1) == i;
%         shapeIx = addInfoVect(:,2) == j;
        
        for k=1:updatedInfo.beam(i).numOfActiveLeafPairs
            % fill the new shapeMap
            vectorIx = updatedInfo.beam(i).shape(j).index + (k-1);
            leftLeafPos = shapeInfoVect(vectorIx);
            rightLeafPos = shapeInfoVect(vectorIx+totalNumOfLeafPairs);
            
            for l=1:size(tempMap,2) % loop over all bixels in this row
                % left bixel edge
                edge_l = updatedInfo.beam(i).posOfCornerBixel(1)...
                        +(l-1-1/2)*updatedInfo.bixelWidth;
                % right bixel edge
                edge_r = updatedInfo.beam(i).posOfCornerBixel(1)...
                        +(l-1+1/2)*updatedInfo.bixelWidth;
                    
                % determine the necessary left/right position for the
                % opening fraction. This will only correspond to the leaf
                % positions if the current bixel is partially covered by
                % the current leaf pair
                
                x_l = max(edge_l,min(leftLeafPos,edge_r));
                x_r = min(edge_r,max(rightLeafPos,edge_l));
                
                % calculate the opening fraction of the bixel
                openingFrac = (x_r-x_l)/updatedInfo.bixelWidth;
                tempMap(k,l) = openingFrac;                
            end
            
%             updatedInfo.beam(i).shape(j).updatedLeft(k) = leftLeafPos;
%             updatedInfo.beam(i).shape(j).updatedRight(k) = rightLeafPos;
            updatedInfo.beam(i).shape(j).leftLeafPos(k) = leftLeafPos;
            updatedInfo.beam(i).shape(j).rightLeafPos(k) = rightLeafPos;
        end
        
        % save the tempMap
        updatedInfo.beam(i).shape(j).shapeMap = tempMap;
        
        %update the shape weight
        updatedInfo.beam(i).shape(j).weight = shapeInfoVect(shapeInd);
        shapeInd = shapeInd +1;
    end
    
end


%% alter unnützer code
% % % vectL = [];
% % % vectR = [];

% loop over all beams
% % % % % for i=1:numel(updatedInfo.beam)
% % % % %     % loop over all shapes
% % % % %     for j=1:updatedInfo.beam(i).numOfShapes
% % % % %         % create a new shapeMap
% % % % %         tempMap = updatedInfo.beam(i).shape(j).shapeMap;
% % % % %         % get beam and shape Ix
% % % % % %         beamIx = addInfoVect(:,1) == i;
% % % % % %         shapeIx = addInfoVect(:,2) == j;
% % % % %         
% % % % %         for k=1:updatedInfo.beam(i).numOfActiveLeafPairs
% % % % % %             leafIx = addInfoVect(:,3) == k;            
% % % % % %             leafPositions = beamIx.*shapeIx.*leafIx.*shapeInfoVect;
% % % % % %             leftIx = find(leaPositions,'first');
% % % % % %             rightIx = find(leafPositions,'last');
% % % % % % %             vectorIx = shapeInfo.beam(i).shape(j).index + (k-1);
% % % % % % %             leftLeafPos = shapeInfoVect(vectorIx);
% % % % % % %             rightLeafPos = shapeInfoVect(vectorIx+totalNumOfLeafPairs);
% % % % % % %             
% % % % % % %             if leftLeafPos>rightLeafPos
% % % % % % %                 warning(['warning! the leaf ' num2str(k) ' of shape ' num2str(j) ' of beam ' num2str(i) ' are overlapping'])
% % % % % % %             end
% % % % % % %             
% % % % % % %             % get shapeMapIx
% % % % % % %             leftLeafInd = round((leftLeafPos-shapeInfo.beam(i).posOfCornerBixel(1))/shapeInfo.bixelWidth)+1;
% % % % % % %             rightLeafInd = round((rightLeafPos-shapeInfo.beam(i).posOfCornerBixel(1))/shapeInfo.bixelWidth)+1;
% % % % % % %             
% % % % % % %             vectL = [vectL; leftLeafInd];
% % % % % % %             vectR = [vectR; rightLeafInd];
% % % % % % %             
% % % % % % %             % set all values in this row to 1
% % % % % % %             tempMap(k,:) = 1;
% % % % % % %             %set all values left of the left leaf edge to zero
% % % % % % %             tempMap(k,1:leftLeafInd-1) = 0;
% % % % % % %             % set all values right of the right leaf edge to zero
% % % % % % %             tempMap(k,rightLeafInd+1:end) = 0;
% % % % % % %             
% % % % % % %             if leftLeafInd == 0
% % % % % % % %                erster bixel offen (zumindest für das linke leaf) 
% %             end

end