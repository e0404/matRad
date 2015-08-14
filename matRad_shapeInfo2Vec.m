function [shapeInfoVec, mappingMx, limMx] = matRad_shapeInfo2Vec(shapeInfo)

% function to create a single vector for the direct aperature optimization
% first: aperature weights
% second: left leaf positions
% third: right leaf positions

% initializing variables

shapeInfoVec = NaN * ones(shapeInfo.totalNumOfShapes+shapeInfo.totalNumOfLeafPairs*2,1);

offset = 0;

%% 1. aperture weights
for i = 1:size(shapeInfo.beam,2)
    for j = 1:shapeInfo.beam(i).numOfShapes
        shapeInfoVec(offset+j) = shapeInfo.beam(i).shape(j).weight;            
    end
    offset = offset + shapeInfo.beam(i).numOfShapes;
end

% 2. left and right leaf positions
%% fill the vector for all shapes of all beams
for i = 1:size(shapeInfo.beam,2)
    for j = 1:shapeInfo.beam(i).numOfShapes
        shapeInfoVec(offset+[1:shapeInfo.beam(i).numOfActiveLeafPairs]) = shapeInfo.beam(i).shape(j).leftLeafPos;
        shapeInfoVec(offset+[1:shapeInfo.beam(i).numOfActiveLeafPairs]+shapeInfo.totalNumOfLeafPairs) = shapeInfo.beam(i).shape(j).rightLeafPos;
        offset = offset + shapeInfo.beam(i).numOfActiveLeafPairs;
    end
end
    

%% 3. create additional information for later use
if nargout > 1
    
    mappingMx =  NaN * ones(shapeInfo.totalNumOfShapes+shapeInfo.totalNumOfLeafPairs*2,3);
    limMx     =  NaN * ones(shapeInfo.totalNumOfShapes+shapeInfo.totalNumOfLeafPairs*2,2);
    limMx(1:shapeInfo.totalNumOfShapes,:) = ones(shapeInfo.totalNumOfShapes,1)*[0 inf];
    
    counter = shapeInfo.totalNumOfShapes + 1;

    for i = 1:numel(shapeInfo.beam)
        for j = 1:shapeInfo.beam(i).numOfShapes
            for k = 1:shapeInfo.beam(i).numOfActiveLeafPairs
                mappingMx(counter,1) = i;
                mappingMx(counter,2) = j;
                mappingMx(counter,3) = k;
                limMx(counter,1)     = shapeInfo.beam(i).lim_l(k);
                limMx(counter,2)     = shapeInfo.beam(i).lim_r(k);
                counter = counter + 1;
            end
        end
    end
    
    mappingMx(counter:end,:) = mappingMx(shapeInfo.totalNumOfShapes+1:counter-1,:);
    limMx(counter:end,:)     = limMx(shapeInfo.totalNumOfShapes+1:counter-1,:);

end
