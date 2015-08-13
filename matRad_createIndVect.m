function indVect = matRad_createIndVect(shapeInfoVect,addInfoVect, shapeInfo)

% function to find the bixel indices corresponding to the values in the
% shapeInfoVect

% initialzing variables
numOfShapes = [shapeInfo.beam.numOfShapes];
totalNumOfShapes = sum(numOfShapes);
totalNumOfLeafPairs = [shapeInfo.beam(:).numOfShapes]*...
                            [shapeInfo.beam(:).numOfActiveLeafPairs]';
indVect = NaN*ones(totalNumOfShapes + totalNumOfLeafPairs,1);

%% 1. Fill the shape weight part
% store the first bixel of the shape
% offset = 0;
% for i=1:numel(shapeInfo.beam)
%     for j=1:numOfShapes(i)
%         indVect(j+offset) = min(shapeInfo.beam(i).bixelIndMap(:));
%     end
%     offset = offset + numOfShapes(i);
% end

% or just leave it as NaN

%% 2. Get the bixel indices from the left and right leaf position

% loop over the remaining part of the shapeInfoVect (leaf positions)
% First all left Leaf positions
for i=totalNumOfShapes+1:totalNumOfShapes+totalNumOfLeafPairs
    % get the information about the the corresponding beam, shape and leaf
    % pair
    currBeam = addInfoVect(i,1);
    currLeaf = addInfoVect(i,3);
    indVect(i) = matRad_getBixelIndForLeaf(shapeInfo,currBeam,currLeaf,shapeInfoVect(i),1);
    if isnan(indVect(i))
        warning(['problem at position ' num2str(i) ', NaN'])
    end
end
% Second all right leaf positions
for i=totalNumOfShapes+totalNumOfLeafPairs+1:(totalNumOfShapes+2*totalNumOfLeafPairs)
    % get the information about the the corresponding beam, shape and leaf
    % pair
    currBeam = addInfoVect(i,1);
    currLeaf = addInfoVect(i,3);
    indVect(i) = matRad_getBixelIndForLeaf(shapeInfo,currBeam,currLeaf,shapeInfoVect(i),0);
    if isnan(indVect(i))
        warning(['problem at position ' num2str(i) ', NaN'])
    end
end
    

end