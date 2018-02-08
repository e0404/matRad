function updatedInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to translate the vector representation of the aperture
% shape and weight into an aperture info struct. At the same time, the
% updated bixel weight vector w is computed and a vector listing the
% correspondence between leaf tips and bixel indices for gradient
% calculation
%
% call
%   [updatedInfo,w,indVect] = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVect)
%
% input
%   apertureInfo:     aperture shape info struct
%   apertureInfoVect: aperture weights and shapes parameterized as vector
%
% output
%   updatedInfo: updated aperture shape info struct according to apertureInfoVect
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

% function to update the apertureInfo struct after the each iteraton of the
% optimization

w = zeros(apertureInfo.totalNumOfBixels,1);

% initializing variables
updatedInfo = apertureInfo;

updatedInfo.apertureVector = apertureInfoVect;

shapeInd = 1;

indVect = NaN*ones(apertureInfo.totalNumOfShapes + apertureInfo.totalNumOfLeafPairs,1);

%% update the shapeMaps
% here the new colimator positions are used to create new shapeMaps that
% now include decimal values instead of binary

% loop over all beams
for i = 1:numel(updatedInfo.beam)
    
    %posOfRightCornerPixel = apertureInfo.beam(i).posOfCornerBixel(1) + (size(apertureInfo.beam(i).bixelIndMap,2)-1)*apertureInfo.bixelWidth;

    % pre compute left and right bixel edges
    edges_l = updatedInfo.beam(i).posOfCornerBixel(1)...
                 + ([1:size(apertureInfo.beam(i).bixelIndMap,2)]-1-1/2)*updatedInfo.bixelWidth;
    edges_r = updatedInfo.beam(i).posOfCornerBixel(1)...
                 + ([1:size(apertureInfo.beam(i).bixelIndMap,2)]-1+1/2)*updatedInfo.bixelWidth;
    
    % loop over all shapes
    for j = 1:updatedInfo.beam(i).numOfShapes
        
        % update the shape weight
        updatedInfo.beam(i).shape(j).weight = apertureInfoVect(shapeInd);
        
        % get dimensions of 2d matrices that store shape/bixel information
        n = apertureInfo.beam(i).numOfActiveLeafPairs;
        
        % extract left and right leaf positions from shape vector
        vectorIx     = updatedInfo.beam(i).shape(j).vectorOffset + ([1:n]-1);
        leftLeafPos  = apertureInfoVect(vectorIx);
        rightLeafPos = apertureInfoVect(vectorIx+apertureInfo.totalNumOfLeafPairs);
        
        % update information in shape structure
        updatedInfo.beam(i).shape(j).leftLeafPos  = leftLeafPos;
        updatedInfo.beam(i).shape(j).rightLeafPos = rightLeafPos;
        
        
        xPosIndLeftLeaf  = round((leftLeafPos - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
        xPosIndRightLeaf = round((rightLeafPos - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
        
        % check limits because of rounding off issues at maximum, i.e., 
        % enfore round(X.5) -> X
        xPosIndLeftLeaf(leftLeafPos == apertureInfo.beam(i).lim_r) = ...
            .5 + (leftLeafPos(leftLeafPos == apertureInfo.beam(i).lim_r) ...
            - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth;
        xPosIndRightLeaf(rightLeafPos == apertureInfo.beam(i).lim_r) = ...
            .5 + (rightLeafPos(rightLeafPos == apertureInfo.beam(i).lim_r) ...
            - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth;

        % find the bixel index that the leaves currently touch
        bixelIndLeftLeaf  = apertureInfo.beam(i).bixelIndMap((xPosIndLeftLeaf-1)*n+[1:n]');
        bixelIndRightLeaf = apertureInfo.beam(i).bixelIndMap((xPosIndRightLeaf-1)*n+[1:n]');
        
        % store information in index vector for gradient calculation
        indVect(apertureInfo.beam(i).shape(j).vectorOffset+[1:n]-1) = bixelIndLeftLeaf;
        indVect(apertureInfo.beam(i).shape(j).vectorOffset+[1:n]-1+apertureInfo.totalNumOfLeafPairs) = bixelIndRightLeaf;

        % calculate opening fraction for every bixel in shape to construct
        % bixel weight vector
        
        coveredByLeftLeaf  = bsxfun(@minus,leftLeafPos,edges_l)  / updatedInfo.bixelWidth;
        coveredByRightLeaf = bsxfun(@minus,edges_r,rightLeafPos) / updatedInfo.bixelWidth;

        tempMap = 1 - (coveredByLeftLeaf  + abs(coveredByLeftLeaf))  / 2 ...
                    - (coveredByRightLeaf + abs(coveredByRightLeaf)) / 2;
            
        % find open bixels
        tempMapIx = tempMap > 0;
        
        currBixelIx = apertureInfo.beam(i).bixelIndMap(tempMapIx);
        w(currBixelIx) = w(currBixelIx) + tempMap(tempMapIx)*updatedInfo.beam(i).shape(j).weight;
    
        % save the tempMap (we need to apply a positivity operator !)
        updatedInfo.beam(i).shape(j).shapeMap = (tempMap  + abs(tempMap))  / 2;
        
        % increment shape index
        shapeInd = shapeInd +1;
    end
    
end

updatedInfo.bixelWeights = w;
updatedInfo.bixelIndices = indVect;

end