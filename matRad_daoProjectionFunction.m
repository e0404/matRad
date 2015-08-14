function [shapeInfoVect, isConstrActive] = matRad_daoProjectionFunction(shapeInfoVect,limVect,totalNumOfShapes)
% function to correct the
% shapeInfoVect[shapeWeights,LeftLeafPos,RightLeafPos] for problems:
% 1. negative weights
% 2. Overlap of left and right leaf
% 3. Leaf out of bounds

% project to feasible set
isConstrActive = shapeInfoVect<=limVect(:,1);
shapeInfoVect(shapeInfoVect<=limVect(:,1)) = limVect(shapeInfoVect<=limVect(:,1),1);

isConstrActive = shapeInfoVect>=limVect(:,2) | isConstrActive;
shapeInfoVect(shapeInfoVect>=limVect(:,2)) = limVect(shapeInfoVect>=limVect(:,2),2);

% correct overlapping leaves
totalNumOfLeafPairs = (numel(shapeInfoVect)-totalNumOfShapes)/2;
leftLeafPos = shapeInfoVect(totalNumOfShapes+[1:totalNumOfLeafPairs]);
rightLeafPos = shapeInfoVect(totalNumOfShapes+totalNumOfLeafPairs+[1:totalNumOfLeafPairs]);

overlapping = leftLeafPos >= rightLeafPos;

leftLeafPos(overlapping)  = mean([leftLeafPos(overlapping) rightLeafPos(overlapping)],2);
rightLeafPos(overlapping) = leftLeafPos(overlapping);


shapeInfoVect(totalNumOfShapes+[1:2*totalNumOfLeafPairs]) = [leftLeafPos;rightLeafPos];


overlapping = [zeros(totalNumOfShapes,1);overlapping;overlapping];
isConstrActive = isConstrActive | overlapping;




end

