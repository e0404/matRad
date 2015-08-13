function [DAOvect, addInfoVect, limVect] = matRad_shapeInfo2Vect(shapeInfo)

% function to create a single vector for the direct aperature optimization
% first: aperature weights
% second: left leaf positions
% third: right leaf positions

% initializing variables
NumOfShapes = [shapeInfo.beam.numOfShapes];
numOfShapes = [shapeInfo.beam.numOfShapes];
totalNumOfShapes = sum(numOfShapes);
totalNumOfLeafPairs = [shapeInfo.beam(:).numOfShapes]*...
                            [shapeInfo.beam(:).numOfActiveLeafPairs]';
apW = zeros(sum(NumOfShapes),1);
%% 1. bixel weights

    
    offset = 0;
    for i = 1:numel(shapeInfo.beam)
        for j = 1:NumOfShapes(i)
            apW(offset+j) = shapeInfo.beam(i).shape(j).weight;            
        end
        offset = offset + NumOfShapes(i);
    end

%% 2. left and right leaf positions
    
    % initializing variables
    leftLeafVect = [];
    rightLeafVect = [];
    
    % fill the vector for all shapes of all beams
    for i=1:numel(shapeInfo.beam)
        for j=1:shapeInfo.beam(i).numOfShapes
            leftLeafVect = [leftLeafVect; shapeInfo.beam(i).shape(j).leftLeafPos];
            rightLeafVect = [rightLeafVect; shapeInfo.beam(i).shape(j).rightLeafPos]; 
        end        
    end    
    
%% 3. create single vector

DAOvect = [apW;leftLeafVect;rightLeafVect];

%% 4. create additional information for later use
    addInfoVect = [];
    
    for i=1:numel(shapeInfo.beam)
        for j=1:shapeInfo.beam(i).numOfShapes
            for k=1:shapeInfo.beam(i).numOfActiveLeafPairs
                currentBeamNum = i;
                currentShapeNum = j; % number within this beam
                currentLeafPairNum = k; % # number within the shape
                addInfoVect = [addInfoVect; currentBeamNum ...
                            currentShapeNum currentLeafPairNum];
            end
        end        
    end
    
addInfoVect = [NaN*ones(size(apW,1),3); addInfoVect; addInfoVect];

%% 5. create limit vector for weight and position limits

limVect = NaN * ones(totalNumOfShapes+2*totalNumOfLeafPairs,2);
% 1. shape weight limits
for i=1:totalNumOfShapes
    limVect(i,:) = [0,inf];
end

% 2. position limits

for i=totalNumOfShapes+1:totalNumOfShapes+2*totalNumOfLeafPairs
    currBeam = addInfoVect(i,1);
    currLeaf = addInfoVect(i,3);
    
    limVect(i,1) = shapeInfo.beam(currBeam).lim_l(currLeaf);
    limVect(i,2) = shapeInfo.beam(currBeam).lim_r(currLeaf);
end

end