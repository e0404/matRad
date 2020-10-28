function [apertureInfoVec, mappingMx, limMx] = matRad_daoApertureInfo2Vec(apertureInfo)
% matRad function to generate a vector respresentation of objects
% A vector representation of the aperture weights and shapes and (optional) 
% some meta information needed during optimization is generated.
%
% call
%   [apertureInfoVec, mappingMx, limMx] = matRad_daoApertureInfo2Vec(apertureInfo)
%
% input
%   apertureInfo:    aperture weight and shape info struct
%
% output
%   apertureInfoVec: vector respresentation of the apertue weights and shapes
%   mappingMx:       mapping of vector components to beams, shapes and leaves
%   limMx:           bounds on vector components, i.e., minimum and maximum
%                    aperture weights (0/inf) and leav positions (custom)
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function to create a single vector for the direct aperature optimization
% first: aperature weights
% second: left leaf positions
% third: right leaf positions

% initializing variables

apertureInfoVec = NaN * ones(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2,1);

offset = 0;

%% 1. aperture weights
for i = 1:size(apertureInfo.beam,2)
    for j = 1:apertureInfo.beam(i).numOfShapes
        apertureInfoVec(offset+j) = apertureInfo.beam(i).shape(j).weight;            
    end
    offset = offset + apertureInfo.beam(i).numOfShapes;
end

% 2. left and right leaf positions
%% fill the vector for all shapes of all beams
for i = 1:size(apertureInfo.beam,2)
    for j = 1:apertureInfo.beam(i).numOfShapes
        apertureInfoVec(offset+[1:apertureInfo.beam(i).numOfActiveLeafPairs]) = apertureInfo.beam(i).shape(j).leftLeafPos;
        apertureInfoVec(offset+[1:apertureInfo.beam(i).numOfActiveLeafPairs]+apertureInfo.totalNumOfLeafPairs) = apertureInfo.beam(i).shape(j).rightLeafPos;
        offset = offset + apertureInfo.beam(i).numOfActiveLeafPairs;
    end
end
    

%% 3. create additional information for later use
if nargout > 1
    
    mappingMx =  NaN * ones(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2,4);
    limMx     =  NaN * ones(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2,2);
    limMx(1:apertureInfo.totalNumOfShapes,:) = ones(apertureInfo.totalNumOfShapes,1)*[0 inf];
    
    counter = 1;
    for i = 1:numel(apertureInfo.beam)
        for j = 1:apertureInfo.beam(i).numOfShapes
            mappingMx(counter,1) = i;
            counter = counter + 1;
        end
    end
    
    shapeOffset = 0;
    for i = 1:numel(apertureInfo.beam)
        for j = 1:apertureInfo.beam(i).numOfShapes
            for k = 1:apertureInfo.beam(i).numOfActiveLeafPairs
                mappingMx(counter,1) = i;
                mappingMx(counter,2) = j + shapeOffset; % store global shape number for grad calc
                mappingMx(counter,3) = j; % store local shape number
                mappingMx(counter,4) = k; % store local leaf number
                
                limMx(counter,1)     = apertureInfo.beam(i).lim_l(k);
                limMx(counter,2)     = apertureInfo.beam(i).lim_r(k);
                counter = counter + 1;
            end
        end
        shapeOffset = shapeOffset + apertureInfo.beam(i).numOfShapes;
    end
    
    mappingMx(counter:end,:) = mappingMx(apertureInfo.totalNumOfShapes+1:counter-1,:);
    limMx(counter:end,:)     = limMx(apertureInfo.totalNumOfShapes+1:counter-1,:);

end
