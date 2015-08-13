function gradVect = matRad_getGradients(wGrad,addInfoVect,indVect,shapeInfo)

% function to determine the Gradients for the leaf positions from the bixel
% weight gradients

% initialize variables
numOfBeams = numel(shapeInfo.beam);
bixelWidth = shapeInfo.bixelWidth;
totalNumOfLeafPairs = [shapeInfo.beam(:).numOfShapes]*...
                            [shapeInfo.beam(:).numOfActiveLeafPairs]';
totalNumOfShapes = sum([shapeInfo.beam(:).numOfShapes]);


aperatureGrad = zeros(totalNumOfShapes,1);
gradVect = NaN * ones(size(addInfoVect,1),1);

%% 1. calculate aperatureGrad
% the gradient of the aperature weight can be calculated from the sum of
% all bixel-gradients of this shape multiplied with the corresponding
% opening fraction

    %  loop over all beams
    offset = 0;
    for i = 1:numOfBeams

        % get the bixelMap
        bixelMap = shapeInfo.beam(i).bixelIndMap;
        firstBixel = min(bixelMap(:));
        lastBixel = max(bixelMap(:));
        
        
        % loop over all shapes
        for j=1:shapeInfo.beam(i).numOfShapes
            
            % add up the gradients x openingFrac for this shape
            for l=firstBixel:lastBixel
                currIx = bixelMap == l;
                aperatureGrad(j+offset) = aperatureGrad(j+offset)+...
                        wGrad(l)*shapeInfo.beam(i).shape(j).shapeMap(currIx);            
            end            
        end
        
        offset = offset + shapeInfo.beam(i).numOfShapes;

    end

%% 2. find the corresponding bixel to the leaf Positions and calculate the gradient

    %loop over all leaf positions
    for i=totalNumOfShapes+1:size(addInfoVect,1)
        % get the beam and shape number for this leaf to find out the
        % aperature weight
        currBeam = addInfoVect(i,1);
        currShape = addInfoVect(i,2);
        
        apWeight = shapeInfo.beam(currBeam).shape(currShape).weight;
        
        % get the bixel gradient from the current Bixel for this leaf
        % position
        currBixel = indVect(i);
        currBixGrad = wGrad(currBixel);
        
        % calculate the gradient for the leaf position
        % g_x = +/- bixelGrad * (aperatureWeight/bixelWidth) the minus sign
        % is used for the left leaf, as a movement in positive x-direction
        % closes the bixel
        gradVect(i) = (currBixGrad * apWeight/bixelWidth);        
        
    end

    % correct the sign for the left leaf positions
    gradVect(totalNumOfShapes+1:totalNumOfShapes+totalNumOfLeafPairs) = ...
        -gradVect(totalNumOfShapes+1:totalNumOfShapes+totalNumOfLeafPairs);

%% 3. form single vector

gradVect(1:length(aperatureGrad)) = aperatureGrad;


end