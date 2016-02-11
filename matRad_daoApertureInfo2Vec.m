function [apertureInfoVec, mappingMx, limMx] = matRad_daoApertureInfo2Vec(apertureInfo)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to generate a vector respresentation of the aperture
% weights and shapes and (optional) some meta information needed during
% optimization
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
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
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
    
    mappingMx =  NaN * ones(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2,2);
    limMx     =  NaN * ones(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2,2);
    limMx(1:apertureInfo.totalNumOfShapes,:) = ones(apertureInfo.totalNumOfShapes,1)*[0 inf];
    
    counter = apertureInfo.totalNumOfShapes + 1;
    shapeOffset = 0;
    for i = 1:numel(apertureInfo.beam)
        for j = 1:apertureInfo.beam(i).numOfShapes
            for k = 1:apertureInfo.beam(i).numOfActiveLeafPairs
                mappingMx(counter,1) = i;
                mappingMx(counter,2) = j + shapeOffset; % store global shape number for grad calc
                %mappingMx(counter,3) = k;
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
