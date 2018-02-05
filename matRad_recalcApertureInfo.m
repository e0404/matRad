function recalc = matRad_recalcApertureInfo(recalc,apertureInfoOld)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to apertures for a different dose resolution
%
% call
%   recalc = matRad_recalcApertureInfo(recalc,apertureInfo)
%
% input
%   recalc:            
%   apertureInfo:
%
% output
%   recalc:            
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

pln = recalc.pln;
stf = recalc.stf;

apertureInfoNew = apertureInfoOld;
apertureInfoNew = rmfield(apertureInfoNew,'beam');

apertureInfoNew.totalNumOfBixels = sum([stf(:).totalNumOfBixels]);

shapeInd = 1;

if recalc.interpNew
    oldGantryAngles = zeros(1,numel(apertureInfoOld.beam));
    oldLeftLeafPoss = zeros(apertureInfoOld.beam(1).numOfActiveLeafPairs,numel(apertureInfoOld.beam));
    oldRightLeafPoss = zeros(apertureInfoOld.beam(1).numOfActiveLeafPairs,numel(apertureInfoOld.beam));
    for i = 1:numel(apertureInfoOld.beam)
        oldGantryAngles(i) = apertureInfoOld.beam(i).gantryAngle;
        oldLeftLeafPoss(:,i) = apertureInfoOld.beam(i).shape(1).leftLeafPos;
        oldRightLeafPoss(:,i) = apertureInfoOld.beam(i).shape(1).rightLeafPos;
    end
end

% MLC parameters:
bixelWidth = stf(1).bixelWidth; % [mm]
numOfMLCLeafPairs = 80;
%     define central leaf pair (here we want the 0mm position to be in the
%     center of a leaf pair (e.g. leaf 41 stretches from -2.5mm to 2.5mm
%     for a bixel/leafWidth of 5mm and 81 leaf pairs)
centralLeafPair = ceil(numOfMLCLeafPairs/2);

% initializing variables
totalNumOfShapes = numel(stf);
for i = 1:numel(apertureInfoOld.beam)
    newInd = (apertureInfoOld.beam(i).doseAngleBorders(1) <= [stf.gantryAngle] & [stf.gantryAngle] <= apertureInfoOld.beam(i).doseAngleBorders(2)).*(1:numel([stf.gantryAngle]));
    newInd(newInd == 0) = [];
    
    for j = newInd
        % get x- and z-coordinates of bixels
        rayPos_bev = reshape([stf(j).ray.rayPos_bev],3,[]);
        X = rayPos_bev(1,:)';
        Z = rayPos_bev(3,:)';
        
        % create ray-map
        maxX = max(X); minX = min(X);
        maxZ = max(Z); minZ = min(Z);
        
        dimX = (maxX-minX)/stf(j).bixelWidth + 1;
        dimZ = (maxZ-minZ)/stf(j).bixelWidth + 1;
        
        rayMap = zeros(dimZ,dimX);
        
        % get indices for x and z positions
        xPos = (X-minX)/stf(j).bixelWidth + 1;
        zPos = (Z-minZ)/stf(j).bixelWidth + 1;
        
        % get indices in the ray-map
        indInRay = zPos + (xPos-1)*dimZ;
        
        % fill ray-map
        rayMap(indInRay) = 1;
        
        % create map of bixel indices
        bixelIndMap = NaN * ones(dimZ,dimX);
        bixelIndMap(indInRay) = [1:stf(j).numOfRays] + (j-1)*stf(1).numOfRays;
        
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
        
        leafPairPos = unique(Z);
        
        % find upmost and downmost leaf pair
        topLeafPairPos = maxZ;
        bottomLeafPairPos = minZ;
        
        topLeafPair = centralLeafPair - topLeafPairPos/bixelWidth;
        bottomLeafPair = centralLeafPair - bottomLeafPairPos/bixelWidth;
        
        % create bool map of active leaf pairs
        isActiveLeafPair = zeros(numOfMLCLeafPairs,1);
        isActiveLeafPair(topLeafPair:bottomLeafPair) = 1;
        
        MLCWindow = [minX-bixelWidth/2 maxX+bixelWidth/2 ...
            minZ-bixelWidth/2 maxZ+bixelWidth/2];
        
        
        % save data for each beam
        apertureInfoNew.beam(j).numOfActiveLeafPairs = dimZ;
        apertureInfoNew.beam(j).leafPairPos = leafPairPos;
        apertureInfoNew.beam(j).isActiveLeafPair = isActiveLeafPair;
        apertureInfoNew.beam(j).centralLeafPair = centralLeafPair;
        apertureInfoNew.beam(j).lim_l = lim_l;
        apertureInfoNew.beam(j).lim_r = lim_r;
        apertureInfoNew.beam(j).bixelIndMap = bixelIndMap;
        apertureInfoNew.beam(j).posOfCornerBixel = posOfCornerBixel;
        apertureInfoNew.beam(j).MLCWindow = MLCWindow;
        apertureInfoNew.beam(j).bixOffset = 1+(j-1)*dimZ;
        apertureInfoNew.beam(j).shape(1).vectorOffset = totalNumOfShapes+1+(j-1)*dimZ;
        
        %inherit from old beam
        apertureInfoNew.beam(j).leafDir = apertureInfoOld.beam(i).leafDir;
        
        %specific to new beam
        apertureInfoNew.beam(j).gantryAngle = pln.gantryAngles(j);
        apertureInfoNew.beam(j).doseAngleBorders = stf(j).doseAngleBorders;
        apertureInfoNew.beam(j).doseAngleBorderCentreDiff = stf(j).doseAngleBorderCentreDiff;
        apertureInfoNew.beam(j).doseAngleBordersDiff = stf(j).doseAngleBordersDiff;
        apertureInfoNew.beam(j).lastOptIndex = stf(j).lastOptIndex;
        apertureInfoNew.beam(j).nextOptIndex = stf(j).lastOptIndex;
        
        
        amountOfOldSpeed = (min(apertureInfoNew.beam(j).doseAngleBorders(2),apertureInfoOld.beam(i).doseAngleBorders(2))-max(apertureInfoNew.beam(j).doseAngleBorders(1),apertureInfoOld.beam(i).doseAngleBorders(1)))./apertureInfoNew.beam(j).doseAngleBordersDiff;
        amountOfOldWeight = (min(apertureInfoNew.beam(j).doseAngleBorders(2),apertureInfoOld.beam(i).doseAngleBorders(2))-max(apertureInfoNew.beam(j).doseAngleBorders(1),apertureInfoOld.beam(i).doseAngleBorders(1)))./apertureInfoOld.beam(i).doseAngleBordersDiff;
        
        amountOfOldWeight_I = (min(apertureInfoNew.beam(j).gantryAngle,apertureInfoOld.beam(i).doseAngleBorders(2))-max(apertureInfoNew.beam(j).doseAngleBorders(1),apertureInfoOld.beam(i).doseAngleBorders(1)))./apertureInfoOld.beam(i).doseAngleBordersDiff;
        amountOfOldWeight_F = (min(apertureInfoNew.beam(j).doseAngleBorders(2),apertureInfoOld.beam(i).doseAngleBorders(2))-max(apertureInfoNew.beam(j).gantryAngle,apertureInfoOld.beam(i).doseAngleBorders(1)))./apertureInfoOld.beam(i).doseAngleBordersDiff;
        
        if ~isfield(apertureInfoNew.beam(j),'gantryRot') || isempty(apertureInfoNew.beam(j).gantryRot)
            apertureInfoNew.beam(j).gantryRot = 0;
            apertureInfoNew.beam(j).shape(1).weight = 0;
            apertureInfoNew.beam(j).shape(1).weight_I = 0;
            apertureInfoNew.beam(j).shape(1).weight_F = 0;
        end
        apertureInfoNew.beam(j).gantryRot = amountOfOldSpeed*apertureInfoOld.beam(i).gantryRot+apertureInfoNew.beam(j).gantryRot;
        
        %recalculate weight, MU
        apertureInfoNew.beam(j).shape(1).weight = apertureInfoNew.beam(j).shape(1).weight+amountOfOldWeight*apertureInfoOld.beam(i).shape(1).weight;
        apertureInfoNew.beam(j).shape(1).weight_I = apertureInfoNew.beam(j).shape(1).weight_I+amountOfOldWeight_I*apertureInfoOld.beam(i).shape(1).weight;
        apertureInfoNew.beam(j).shape(1).weight_F = apertureInfoNew.beam(j).shape(1).weight_F+amountOfOldWeight_F*apertureInfoOld.beam(i).shape(1).weight;
        apertureInfoNew.beam(j).MU = apertureInfoNew.beam(j).shape(1).weight.*apertureInfoNew.weightToMU;
        
        apertureInfoNew.beam(j).MURate = apertureInfoNew.beam(j).MU.*apertureInfoNew.beam(j).gantryRot./apertureInfoNew.beam(j).doseAngleBordersDiff;
        
        %apertureInfoNew.beam(j).shape(1).jacobiScale = apertureInfoOld.beam(i).shape(1).jacobiScale;
        apertureInfoNew.beam(j).shape(1).jacobiScale = 1;
        
        if recalc.interpNew
            %interpolate new apertures now so that weights are not
            %overwritten
            apertureInfoNew.beam(j).shape(1).leftLeafPos = (interp1(oldGantryAngles',oldLeftLeafPoss',apertureInfoNew.beam(j).gantryAngle))';
            apertureInfoNew.beam(j).shape(1).rightLeafPos = (interp1(oldGantryAngles',oldRightLeafPoss',apertureInfoNew.beam(j).gantryAngle))';
            
            apertureInfoNew.beam(j).shape(1).leftLeafPos_I = (interp1(oldGantryAngles',oldLeftLeafPoss',apertureInfoNew.beam(j).doseAngleBorders(1)))';
            apertureInfoNew.beam(j).shape(1).rightLeafPos_I = (interp1(oldGantryAngles',oldRightLeafPoss',apertureInfoNew.beam(j).doseAngleBorders(1)))';
            
            apertureInfoNew.beam(j).shape(1).leftLeafPos_F = (interp1(oldGantryAngles',oldLeftLeafPoss',apertureInfoNew.beam(j).doseAngleBorders(2)))';
            apertureInfoNew.beam(j).shape(1).rightLeafPos_F = (interp1(oldGantryAngles',oldRightLeafPoss',apertureInfoNew.beam(j).doseAngleBorders(2)))';
        else
            apertureInfoNew.beam(j).shape(1).leftLeafPos = apertureInfoOld.beam(i).shape(1).leftLeafPos;
            apertureInfoNew.beam(j).shape(1).rightLeafPos = apertureInfoOld.beam(i).shape(1).rightLeafPos;
        end
        
        %all beams are now "optimized" to prevent their weights from being
        %overwritten
        %optAngleBorders becomes doseAngleBorders
        apertureInfoNew.beam(j).numOfShapes = 1;
        apertureInfoNew.beam(j).optimizeBeam = true;
        apertureInfoNew.beam(j).doseAngleOpt = stf(j).doseAngleOpt;
        apertureInfoNew.beam(j).optAngleBorders = stf(j).doseAngleBorders;
        apertureInfoNew.beam(j).optAngleBorderCentreDiff = stf(j).doseAngleBorderCentreDiff;
        apertureInfoNew.beam(j).optAngleBordersDiff = stf(j).doseAngleBordersDiff;
        apertureInfoNew.beam(j).timeFacCurr = stf(j).timeFacCurr;
        apertureInfoNew.beam(j).timeFacPrev = stf(j).timeFacPrev;
        apertureInfoNew.beam(j).timeFacNext = stf(j).timeFacNext;
        apertureInfoNew.beam(j).IandFTimeInd = stf(j).IandFTimeInd;
        apertureInfoNew.beam(j).IandFFac = stf(j).IandFFac;
        apertureInfoNew.beam(j).timeFac = stf(j).timeFac;
        
        %{
        if isfield(recalc,'dijNew') && ~recalc.dijNew
            oldGantryAngles = [apertureInfoOld.beam.gantryAngle];
            diff = abs(apertureInfoNew.beam(j).gantryAngle-[apertureInfoOld.beam.gantryAngle]);
            minDiffInd = diff == min(diff);
            if isempty(stf(j).copyInd)
                recalc.pln.gantryAngles(j) = oldGantryAngles(minDiffInd);
            elseif stf(j).copyInd == 1
                recalc.pln.gantryAngles(j) = min(oldGantryAngles(minDiffInd));
            elseif stf(j).copyInd == 2
                recalc.pln.gantryAngles(j) = max(oldGantryAngles(minDiffInd));
            end
        end
        %}
        
        if apertureInfoOld.beam(i).initializeBeam
            apertureInfoNew.beam(j).initializeBeam = true;
            apertureInfoNew.beam(j).initAngleBorders = stf(j).initAngleBorders;
            apertureInfoNew.beam(j).initAngleBorderCentreDiff = stf(j).initAngleBorderCentreDiff;
            apertureInfoNew.beam(j).initAngleBordersDiff = stf(j).initAngleBordersDiff;
        else
            apertureInfoNew.beam(j).initializeBeam = false;
        end
        
        apertureInfoNew.apertureVector(shapeInd) = apertureInfoNew.beam(j).shape(1).weight;
        shapeInd = shapeInd+1;
    end
end

apertureInfoNew.totalNumOfShapes = sum([apertureInfoNew.beam.numOfShapes]);
apertureInfoNew.totalNumOfLeafPairs = sum([apertureInfoNew.beam.numOfShapes]*[apertureInfoNew.beam.numOfActiveLeafPairs]');
apertureInfoNew.doseTotalNumOfLeafPairs = sum([apertureInfoNew.beam(:).numOfActiveLeafPairs]);

%recalc apertureVector
[apertureInfoNew.apertureVector, apertureInfoNew.mappingMx, apertureInfoNew.limMx] = matRad_daoApertureInfo2Vec(apertureInfoNew);

recalc.apertureInfo = apertureInfoNew;
recalc.stf = stf;