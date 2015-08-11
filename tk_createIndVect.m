function indVect = tk_createIndVect(shapeInfoVect,addInfoVect, shapeInfo)

% function to find the bixel indices corresponding to the values in the
% shapeInfoVect

% initialzing variables
indVect = NaN*ones(shapeInfo.totalNumOfBixels,1);
totalNumOfLeafPairs = [shapeInfo.beam(:).numOfShapes]*...
                            [shapeInfo.beam(:).numOfActiveLeafPairs]';

%% 1. Fill the bixel weight part
indVect(1:shapeInfo.totalNumOfBixels) = 1:shapeInfo.totalNumOfBixels;

%% 2. Get the bixel indices from the left and right leaf position

% loop over the remaining part of the shapeInfoVect (leaf positions)
% First all left Leaf positions
for i=shapeInfo.totalNumOfBixels+1:numel(shapeInfoVect)-totalNumOfLeafPairs
    % get the information about the the corresponding beam, shape and leaf
    % pair
    currBeam = addInfoVect(i,1);
    currLeaf = addInfoVect(i,3);
    indVect(i) = tk_getBixelInd(shapeInfo,currBeam,currLeaf,shapeInfoVect(i),1);
    if isnan(indVect(i))
        warning(['problem at position ' num2str(i) ', NaN'])
    end
end
% Second all right leaf positions
for i=numel(shapeInfoVect)-totalNumOfLeafPairs+1:numel(shapeInfoVect)
    % get the information about the the corresponding beam, shape and leaf
    % pair
    currBeam = addInfoVect(i,1);
    currLeaf = addInfoVect(i,3);
    indVect(i) = tk_getBixelInd(shapeInfo,currBeam,currLeaf,shapeInfoVect(i),0);
    if isnan(indVect(i))
        warning(['problem at position ' num2str(i) ', NaN'])
    end
end
    

end