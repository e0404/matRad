function apertureInfo = matRad_sequencing2ApertureInfo(sequencing,stf,pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to generate a shape info struct based on the result of
% multileaf collimator sequencing
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
totalNumOfShapes = sum([sequencing.beam.numOfShapes]);
weightOffset = 1;
vectorOffset = totalNumOfShapes + 1; % used for bookkeeping in the vector for optimization
bixOffset = 1; %used for gradient calculations
interpGetsTransition = false; % boolean to determine if an interpolated beam is responsible for a leaf speed constraint check

if pln.propOpt.runVMAT
    totalNumOfOptBixels = 0;
    totalNumOfLeafPairs = 0;
    
    apertureInfo.propVMAT.jacobT = zeros(sum([sequencing.beam.numOfShapes]),numel(sequencing.beam));
end

apertureInfo.jacobiScale = ones(totalNumOfShapes,1);

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
    bixelIndMap(indInRay) = (1:stf(i).numOfRays) + bixelIndOffset;
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
    for m = 1:sequencing.beam(i).numOfShapes
        
        if isfield(sequencing.beam(i),'shapes')
            
            % loading shape from Sequencing result
            shapeMap = sequencing.beam(i).shapes(:,:,m);
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
                    
                    %Can happen in some cases in SW trajectory sampling
                    if leftLeafPos(l) < lim_l(l)
                        leftLeafPos(l) = lim_l(l);
                    end
                    if rightLeafPos(l) > lim_r(l)
                        rightLeafPos(l) = lim_r(l);
                    end
                    
                end
            end
            
            % save data for each shape of this beam
            apertureInfo.beam(i).shape(m).leftLeafPos = leftLeafPos;
            apertureInfo.beam(i).shape(m).rightLeafPos = rightLeafPos;
            apertureInfo.beam(i).shape(m).weight = sequencing.beam(i).shapesWeight(m);
            apertureInfo.beam(i).shape(m).shapeMap = shapeMap;
        end
        
        if pln.propOpt.runVMAT
            apertureInfo.beam(i).shape(m).MURate = sequencing.beam(i).MURate;
        end
        
        apertureInfo.beam(i).shape(m).jacobiScale = 1;
        
        if pln.propOpt.runVMAT && pln.propOpt.VMAToptions.continuousAperture
            apertureInfo.beam(i).shape(m).vectorOffset = [vectorOffset vectorOffset+dimZ];
            
            % update index for bookkeeping
            vectorOffset = vectorOffset + dimZ*nnz(stf(i).propVMAT.doseAngleDAO);
        else
            apertureInfo.beam(i).shape(m).vectorOffset = vectorOffset;
            
            % update index for bookkeeping
            vectorOffset = vectorOffset + dimZ;
        end
        apertureInfo.beam(i).shape(m).weightOffset = weightOffset;
        weightOffset = weightOffset+1;
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
    apertureInfo.beam(i).numOfShapes = sequencing.beam(i).numOfShapes;
    apertureInfo.beam(i).numOfActiveLeafPairs = dimZ;
    apertureInfo.beam(i).leafPairPos = leafPairPos;
    apertureInfo.beam(i).isActiveLeafPair = isActiveLeafPair;
    apertureInfo.beam(i).centralLeafPair = centralLeafPair;
    apertureInfo.beam(i).lim_l = lim_l;
    apertureInfo.beam(i).lim_r = lim_r;
    apertureInfo.beam(i).bixelIndMap = bixelIndMap;
    apertureInfo.beam(i).posOfCornerBixel = posOfCornerBixel;
    apertureInfo.beam(i).MLCWindow = MLCWindow;
    apertureInfo.beam(i).gantryAngle = stf(i).gantryAngle;
    
    if pln.propOpt.runVMAT
        
        apertureInfo.beam(i).bixOffset = bixOffset;
        bixOffset = bixOffset+apertureInfo.beam(i).numOfActiveLeafPairs;
        
        apertureInfo.propVMAT.beam(i).DAOBeam = stf(i).propVMAT.DAOBeam;
        apertureInfo.propVMAT.beam(i).FMOBeam = stf(i).propVMAT.FMOBeam;
        
        apertureInfo.propVMAT.beam(i).leafDir = sequencing.beam(i).leafDir;
        
        apertureInfo.propVMAT.beam(i).doseAngleBorders = stf(i).propVMAT.doseAngleBorders;
        apertureInfo.propVMAT.beam(i).doseAngleBorderCentreDiff = stf(i).propVMAT.doseAngleBorderCentreDiff;
        apertureInfo.propVMAT.beam(i).doseAngleBordersDiff = stf(i).propVMAT.doseAngleBordersDiff;
        
        if apertureInfo.propVMAT.beam(i).DAOBeam
            
            totalNumOfOptBixels = totalNumOfOptBixels+stf(i).totalNumOfBixels;
            totalNumOfLeafPairs = totalNumOfLeafPairs+apertureInfo.beam(i).numOfShapes*apertureInfo.beam(i).numOfActiveLeafPairs;
            
            apertureInfo.beam(i).gantryRot = sequencing.beam(i).gantryRot;
            
            apertureInfo.propVMAT.beam(i).DAOAngleBorders = stf(i).propVMAT.DAOAngleBorders;
            apertureInfo.propVMAT.beam(i).DAOAngleBorderCentreDiff = stf(i).propVMAT.DAOAngleBorderCentreDiff;
            apertureInfo.propVMAT.beam(i).DAOAngleBordersDiff = stf(i).propVMAT.DAOAngleBordersDiff;
            
            apertureInfo.propVMAT.beam(i).timeFacCurr = stf(i).propVMAT.timeFacCurr;
            apertureInfo.propVMAT.beam(i).timeFac = stf(i).propVMAT.timeFac;
            
            apertureInfo.propVMAT.beam(i).lastDAOIndex = stf(i).propVMAT.lastDAOIndex;
            apertureInfo.propVMAT.beam(i).nextDAOIndex = stf(i).propVMAT.nextDAOIndex;
            apertureInfo.propVMAT.beam(i).DAOIndex = stf(i).propVMAT.DAOIndex;
            
            if apertureInfo.propVMAT.beam(i).FMOBeam
                apertureInfo.propVMAT.beam(i).FMOAngleBorders = stf(i).propVMAT.FMOAngleBorders;
                apertureInfo.propVMAT.beam(i).FMOAngleBorderCentreDiff = stf(i).propVMAT.FMOAngleBorderCentreDiff;
                apertureInfo.propVMAT.beam(i).FMOAngleBordersDiff = stf(i).propVMAT.FMOAngleBordersDiff;
            end
            
            if pln.propOpt.VMAToptions.continuousAperture
                apertureInfo.propVMAT.beam(i).timeFacInd = stf(i).propVMAT.timeFacInd;
                apertureInfo.propVMAT.beam(i).doseAngleDAO = stf(i).propVMAT.doseAngleDAO;
                
                apertureInfo.propVMAT.beam(i).leafConstMask = 1;
                interpGetsTransition = apertureInfo.propVMAT.beam(i).timeFac(3) ~= 0;
            end
            
            apertureInfo.propVMAT.jacobT(stf(i).propVMAT.DAOIndex,i) = stf(i).propVMAT.timeFacCurr;
            
        else
            apertureInfo.propVMAT.beam(i).fracFromLastDAO = stf(i).propVMAT.fracFromLastDAO;
            apertureInfo.propVMAT.beam(i).timeFracFromLastDAO = stf(i).propVMAT.timeFracFromLastDAO;
            apertureInfo.propVMAT.beam(i).timeFracFromNextDAO = stf(i).propVMAT.timeFracFromNextDAO;
            apertureInfo.propVMAT.beam(i).lastDAOIndex = stf(i).propVMAT.lastDAOIndex;
            apertureInfo.propVMAT.beam(i).nextDAOIndex = stf(i).propVMAT.nextDAOIndex;
            
            if pln.propOpt.VMAToptions.continuousAperture
                apertureInfo.propVMAT.beam(i).fracFromLastDAO_I = stf(i).propVMAT.fracFromLastDAO_I;
                apertureInfo.propVMAT.beam(i).fracFromLastDAO_F = stf(i).propVMAT.fracFromLastDAO_F;
                apertureInfo.propVMAT.beam(i).fracFromNextDAO_I = stf(i).propVMAT.fracFromNextDAO_I;
                apertureInfo.propVMAT.beam(i).fracFromNextDAO_F = stf(i).propVMAT.fracFromNextDAO_F;
            end
            
            apertureInfo.propVMAT.jacobT(stf(stf(i).propVMAT.lastDAOIndex).propVMAT.DAOIndex,i) = stf(stf(i).propVMAT.lastDAOIndex).propVMAT.timeFacCurr.*stf(i).propVMAT.timeFracFromLastDAO.*stf(i).propVMAT.doseAngleBordersDiff./stf(stf(i).propVMAT.lastDAOIndex).propVMAT.doseAngleBordersDiff;
            apertureInfo.propVMAT.jacobT(stf(stf(i).propVMAT.nextDAOIndex).propVMAT.DAOIndex,i) = stf(stf(i).propVMAT.nextDAOIndex).propVMAT.timeFacCurr.*stf(i).propVMAT.timeFracFromNextDAO.*stf(i).propVMAT.doseAngleBordersDiff./stf(stf(i).propVMAT.lastDAOIndex).propVMAT.doseAngleBordersDiff;
            
            if interpGetsTransition
                apertureInfo.propVMAT.beam(i).leafConstMask = 1;
            end
            interpGetsTransition = false;
        end
    end
end

% save global data
apertureInfo.runVMAT            = pln.propOpt.runVMAT;
apertureInfo.preconditioner     = pln.propOpt.preconditioner;
apertureInfo.bixelWidth         = bixelWidth;
apertureInfo.numOfMLCLeafPairs  = numOfMLCLeafPairs;
apertureInfo.totalNumOfBixels   = totalNumOfBixels;
apertureInfo.totalNumOfShapes   = sum([apertureInfo.beam.numOfShapes]);

if isfield(sequencing,'weightToMU')
    apertureInfo.weightToMU = sequencing.weightToMU;
end

if pln.propOpt.runVMAT
    
    tempStruct_beam = apertureInfo.propVMAT.beam;
    tempStruct_jacobT = apertureInfo.propVMAT.jacobT;
    apertureInfo.propVMAT = pln.propOpt.VMAToptions;
    apertureInfo.propVMAT.beam = tempStruct_beam;
    apertureInfo.propVMAT.jacobT = tempStruct_jacobT;
    
    apertureInfo.totalNumOfOptBixels = totalNumOfOptBixels;
    apertureInfo.doseTotalNumOfLeafPairs = sum([apertureInfo.beam(:).numOfActiveLeafPairs]);
    
    if apertureInfo.propVMAT.continuousAperture
        apertureInfo.totalNumOfLeafPairs = sum(reshape([apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).doseAngleDAO],2,[]),1)*[apertureInfo.beam([apertureInfo.propVMAT.beam.DAOBeam]).numOfActiveLeafPairs]';
        
        % count number of transitions
        apertureInfo.propVMAT.numLeafSpeedConstraint      = nnz([apertureInfo.propVMAT.beam.leafConstMask]);
        apertureInfo.propVMAT.numLeafSpeedConstraintDAO   = nnz([apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).leafConstMask]);
    else
        apertureInfo.totalNumOfLeafPairs = totalNumOfLeafPairs;
    end
    
    % fix instances of leaf touching
    apertureInfo = matRad_leafTouching(apertureInfo);
    
    shapeInd = 0;
    for i = 1:numel(apertureInfo.beam)
        if apertureInfo.propVMAT.beam(i).DAOBeam
            shapeInd = shapeInd+1;
            apertureInfo.propVMAT.beam(i).timeInd = apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+shapeInd;
        end
    end
    
else
    apertureInfo.totalNumOfOptBixels = sum(stf(i).totalNumOfBixels);
    apertureInfo.totalNumOfLeafPairs = sum([apertureInfo.beam.numOfShapes]*[apertureInfo.beam.numOfActiveLeafPairs]');
    apertureInfo.doseTotalNumOfLeafPairs = apertureInfo.totalNumOfLeafPairs;
end

% create vectors for optimization
[apertureInfo.apertureVector, apertureInfo.mappingMx, apertureInfo.limMx] = matRad_daoApertureInfo2Vec(apertureInfo);

end






