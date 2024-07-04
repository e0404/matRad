function apertureInfo = matRad_sequencing2ApertureInfo(Sequencing,stf)
% matRad function to generate a shape info struct 
% based on the result of multileaf collimator sequencing
%
% call
%   apertureInfo = matRad_sequencing2ApertureInfo(Sequencing,stf)
%
% input
%   Sequencing: matRad sequencing result struct
%   stf:        matRad steering information struct
%
% output
%   apertureInfo: matRad aperture weight and shape info struct
%
% References
%   
%   -
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

% MLC parameters:
bixelWidth = stf(1).bixelWidth; % [mm]
numOfMLCLeafPairs = 80;
%     define central leaf pair (here we want the 0mm position to be in the
%     center of a leaf pair (e.g. leaf 41 stretches from -2.5mm to 2.5mm
%     for a bixel/leafWidth of 5mm and 81 leaf pairs)
centralLeafPair = ceil(numOfMLCLeafPairs/2);

% initializing variables
bixelIndOffset = 0; % used for creation of bixel index maps
totalNumOfBixels = sum([stf(:).totalNumOfBixels]);
totalNumOfShapes = sum([Sequencing.beam.numOfShapes]);
vectorOffset = totalNumOfShapes + 1; % used for bookkeeping in the vector for optimization

% loop over all beams
for i=1:size(stf,2)
    
    %% 1. read stf and create maps (Ray & Bixelindex)
    
    % get x- and z-coordinates of bixels
    rayPos_bev = reshape([stf(i).ray.rayPos_bev],3,[]);
    X = rayPos_bev(1,:)';
    Z = rayPos_bev(3,:)';
    
    % create ray-map
    maxX = max(X); minX = min(X);    
    maxZ = max(Z); minZ = min(Z);
    
    dimX = (maxX-minX)/stf(i).bixelWidth + 1;
    dimZ = (maxZ-minZ)/stf(i).bixelWidth + 1;

    rayMap = zeros(dimZ,dimX);
    
    % get indices for x and z positions
    xPos = (X-minX)/stf(i).bixelWidth + 1;
    zPos = (Z-minZ)/stf(i).bixelWidth + 1;
    
    % get indices in the ray-map
    indInRay = zPos + (xPos-1)*dimZ;

    % fill ray-map
    rayMap(indInRay) = 1;
    
    % create map of bixel indices
    bixelIndMap = NaN * ones(dimZ,dimX);
    bixelIndMap(indInRay) = [1:stf(i).numOfRays] + bixelIndOffset;
    bixelIndOffset = bixelIndOffset + stf(i).numOfRays;
    
    % store physical position of first entry in bixelIndMap
    posOfCornerBixel = [minX 0 minZ];
    
    % get leaf limits from the leaf map
    lim_l = NaN * ones(dimZ,1);
    lim_r = NaN * ones(dimZ,1);
    % looping oder leaf pairs
    for l = 1:dimZ
        lim_lInd = find(rayMap(l,:),1,'first');
        lim_rInd = find(rayMap(l,:),1,'last');
        % the physical position [mm] can be calculated from the indices
        lim_l(l) = (lim_lInd-1)*bixelWidth + minX - 1/2*bixelWidth;
        lim_r(l) = (lim_rInd-1)*bixelWidth + minX + 1/2*bixelWidth;
    end
        
    % get leaf positions for all shapes
    % leaf positions can be extracted from the shapes created in Sequencing
    for m = 1:Sequencing.beam(i).numOfShapes
        
        % loading shape from Sequencing result
        shapeMap = Sequencing.beam(i).shapes(:,:,m);
        % get left and right leaf indices from shapemap
        % initializing limits
        leftLeafPos = NaN * ones(dimZ,1);
        rightLeafPos = NaN * ones(dimZ,1);
        % looping over leaf pairs
        for l = 1:dimZ
            leftLeafPosInd  = find(shapeMap(l,:),1,'first');
            rightLeafPosInd = find(shapeMap(l,:),1,'last');
            
            if isempty(leftLeafPosInd) && isempty(rightLeafPosInd) % if no bixel is open, use limits from Ray positions
                leftLeafPos(l) = (lim_l(l)+lim_r(l))/2;
                rightLeafPos(l) = leftLeafPos(l);
            else
            % the physical position [mm] can be calculated from the indices
                leftLeafPos(l) = (leftLeafPosInd-1)*bixelWidth...
                                    + minX - 1/2*bixelWidth;
                rightLeafPos(l) = (rightLeafPosInd-1)*bixelWidth...
                                    + minX + 1/2*bixelWidth;
                              
            end
        end
        
        % save data for each shape of this beam
        apertureInfo.beam(i).shape(m).leftLeafPos = leftLeafPos;
        apertureInfo.beam(i).shape(m).rightLeafPos = rightLeafPos;
        apertureInfo.beam(i).shape(m).weight = Sequencing.beam(i).shapesWeight(m);
        apertureInfo.beam(i).shape(m).shapeMap = shapeMap;
        apertureInfo.beam(i).shape(m).vectorOffset = vectorOffset;
        
        % update index for bookkeeping
        vectorOffset = vectorOffset + dimZ;
           
    end
        
    % z-coordinates of active leaf pairs        
    % get z-coordinates from bixel positions
    leafPairPos = unique(Z);
        
    % find upmost and downmost leaf pair
    topLeafPairPos = maxZ;
    bottomLeafPairPos = minZ;
    
    topLeafPair = centralLeafPair - topLeafPairPos/bixelWidth;
    bottomLeafPair = centralLeafPair - bottomLeafPairPos/bixelWidth;
        
    % create bool map of active leaf pairs
    isActiveLeafPair = zeros(numOfMLCLeafPairs,1);
    isActiveLeafPair(topLeafPair:bottomLeafPair) = 1;
        
    % create MLC window
    % getting the dimensions of the MLC in order to be able to plot the
    % shapes using physical coordinates
    MLCWindow = [minX-bixelWidth/2 maxX+bixelWidth/2 ...
                    minZ-bixelWidth/2 maxZ+bixelWidth/2];
    
    % save data for each beam
    apertureInfo.beam(i).numOfShapes = Sequencing.beam(i).numOfShapes;
    apertureInfo.beam(i).numOfActiveLeafPairs = dimZ;
    apertureInfo.beam(i).leafPairPos = leafPairPos;
    apertureInfo.beam(i).isActiveLeafPair = isActiveLeafPair;
    apertureInfo.beam(i).centralLeafPair = centralLeafPair;
    apertureInfo.beam(i).lim_l = lim_l;
    apertureInfo.beam(i).lim_r = lim_r;
    apertureInfo.beam(i).bixelIndMap = bixelIndMap;
    apertureInfo.beam(i).posOfCornerBixel = posOfCornerBixel;
    apertureInfo.beam(i).MLCWindow = MLCWindow;
    
end

% save global data
apertureInfo.bixelWidth = bixelWidth;
apertureInfo.numOfMLCLeafPairs = numOfMLCLeafPairs;
apertureInfo.totalNumOfBixels = totalNumOfBixels;
apertureInfo.totalNumOfShapes = sum([apertureInfo.beam.numOfShapes]);
apertureInfo.totalNumOfLeafPairs = sum([apertureInfo.beam.numOfShapes]*[apertureInfo.beam.numOfActiveLeafPairs]');

% create vectors for optimization
[apertureInfo.apertureVector, apertureInfo.mappingMx, apertureInfo.limMx] = matRad_OptimizationProblemDAO.matRad_daoApertureInfo2Vec(apertureInfo);

end
