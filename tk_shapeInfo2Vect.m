function [DAOvect, addInfoVect] = tk_shapeInfo2Vect(shapeInfo,dij)

% function to create a single vector for the direct aperature optimization
% first: bixel weights
% second: left leaf positions
% third: right leaf positions

%% 1. bixel weights

    % initializing variables
    totalNumOfBixels = shapeInfo.totalNumOfBixels;
    beamNumVect = dij.beamNum;
    w = zeros(totalNumOfBixels,1);
    
    for i = 1:totalNumOfBixels
        % find the bixelposition corresponding to the bixel i

            % 1. find the corresponding beam
            beamNum = beamNumVect(i);

            % 2. find index in the MLC map
            MLCPosInd = shapeInfo.beam(beamNum).bixelIndMap == i;

        % add up the fluence from every shape of this beam
            for j=1:shapeInfo.beam(beamNum).numOfShapes
                w(i) = w(i) + shapeInfo.beam(beamNum).shape(j).weight * ...
                        shapeInfo.beam(beamNum).shape(j).shapeMap(MLCPosInd);
            end
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

DAOvect = [w;leftLeafVect;rightLeafVect];

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
    
addInfoVect = [NaN*ones(size(w,1),3); addInfoVect; addInfoVect];

end