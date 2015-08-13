function [DAOvect, addInfoVect] = tk_shapeInfo2Vect(shapeInfo)

% function to create a single vector for the direct aperature optimization
% first: aperature weights
% second: left leaf positions
% third: right leaf positions

%% 1. bixel weights

    % initializing variables
    NumOfShapes = [shapeInfo.beam.numOfShapes];
    apW = zeros(sum(NumOfShapes),1);
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

end