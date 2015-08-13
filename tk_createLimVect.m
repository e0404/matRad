function limVect = tk_createLimVect(shapeInfo,addInfoVect)

% function to create a vector/matrix containing the limits for:
% 1. shape weigth (>=0, <inf)
% 2. left leafs (>= lim_l, <= lim_r)
% 3. right leafs (>= lim_l, <= lim_r)


% initializing variables
numOfShapes = [shapeInfo.beam.numOfShapes];
totalNumOfShapes = sum(numOfShapes);
totalNumOfLeafPairs = [shapeInfo.beam(:).numOfShapes]*...
                            [shapeInfo.beam(:).numOfActiveLeafPairs]';
limVect = NaN * ones(totalNumOfShapes+2*totalNumOfLeafPairs,2);
%% 1. shape weight limits
for i=1:totalNumOfShapes
    limVect(i,:) = [0,inf];
end

%% 2. limits

for i=totalNumOfShapes+1:totalNumOfShapes+2*totalNumOfLeafPairs
    currBeam = addInfoVect(i,1);
    currLeaf = addInfoVect(i,3);
    
    limVect(i,1) = shapeInfo.beam(currBeam).lim_l(currLeaf);
    limVect(i,2) = shapeInfo.beam(currBeam).lim_r(currLeaf);
end



end