function bixelInd = tk_getBixelInd(shapeInfo,beamNum,leafNum,xPos)

% function to find the bixel index from the information of the leaf
% position of a certain beam
% shapeInfo: struct containing info from sequencing
% beamNum: beam number
% leafNum: number of the leaf of this beam, here the enumeration is from 1
% to numberOfLeaves for this beam. Attention leafNum=1 corresponds to
% lowest leafPair when regarding the physical position

%% 1. store data for chosen beam
try
    info = shapeInfo.beam(beamNum);
catch
    error('invalid beam number or shape structure. Please check the variables!')
end
%% 2. check if the info is valid

if leafNum < 1 || leafNum > info.numOfActiveLeafPairs
    error('specified leaf is not active... please check the value')
end

if xPos < info.lim_l(leafNum) || xPos > info.lim_r(leafNum)
    error('Out of bounds. The specified position is not covered by a ray!')
end

%% 3. find position within the map of the MLC

zPos = info.leafPairPos(leafNum);

zPosInd = (zPos-info.posOfCornerBixel(3))/shapeInfo.bixelWidth+1;
xPosInd = round((xPos-info.posOfCornerBixel(1))/shapeInfo.bixelWidth)+1;

%% 4. get bixel index

bixelInd = info.bixelIndMap(zPosInd,xPosInd);

end