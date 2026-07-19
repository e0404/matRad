function [geometry, bixelIndOffset] = matRad_getMLCGeometry(stfBeam, numOfMLCLeafPairs, bixelIndOffset)
% matRad function to derive the per-beam MLC geometry from an stf beam's
% ray positions: bixel index map, leaf position limits, active leaf pairs
% and MLC window. Shared by aperture-info construction in the photon
% sequencers and the fine-angle recalculation.
%
% call:
%   [geometry, bixelIndOffset] = matRad_getMLCGeometry(stfBeam, numOfMLCLeafPairs, bixelIndOffset)
%
% input:
%   stfBeam:            single stf beam entry (stf(i))
%   numOfMLCLeafPairs:  total number of physical MLC leaf pairs
%   bixelIndOffset:     index of the last bixel of the preceding beams;
%                       the beam's rays are numbered starting from
%                       bixelIndOffset + 1
%
% output:
%   geometry:        struct with fields numOfActiveLeafPairs, leafPairPos,
%                    isActiveLeafPair, centralLeafPair, lim_l, lim_r,
%                    bixelIndMap, posOfCornerBixel, MLCWindow
%   bixelIndOffset:  updated offset (input + stfBeam.numOfRays)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bixelWidth = stfBeam.bixelWidth; % [mm]

% define central leaf pair (here we want the 0mm position to be in the
% center of a leaf pair, e.g. leaf 41 stretches from -2.5mm to 2.5mm
% for a bixel/leafWidth of 5mm and 81 leaf pairs)
centralLeafPair = ceil(numOfMLCLeafPairs / 2);

% get x- and z-coordinates of bixels
rayPos_bev = reshape([stfBeam.ray.rayPos_bev], 3, []);
X = rayPos_bev(1, :)';
Z = rayPos_bev(3, :)';

% create ray-map
maxX = max(X);
minX = min(X);
maxZ = max(Z);
minZ = min(Z);

dimX = (maxX - minX) / bixelWidth + 1;
dimZ = (maxZ - minZ) / bixelWidth + 1;

rayMap = zeros(dimZ, dimX);

% get indices for x and z positions
xPos = (X - minX) / bixelWidth + 1;
zPos = (Z - minZ) / bixelWidth + 1;

% get indices in the ray-map
indInRay = zPos + (xPos - 1) * dimZ;

% fill ray-map
rayMap(indInRay) = 1;

% create map of bixel indices
bixelIndMap = NaN * ones(dimZ, dimX);
bixelIndMap(indInRay) = (1:stfBeam.numOfRays) + bixelIndOffset;
bixelIndOffset = bixelIndOffset + stfBeam.numOfRays;

% get leaf limits from the leaf map
lim_l = NaN * ones(dimZ, 1);
lim_r = NaN * ones(dimZ, 1);
for l = 1:dimZ
    lim_lInd = find(rayMap(l, :), 1, 'first');
    lim_rInd = find(rayMap(l, :), 1, 'last');
    % the physical position [mm] can be calculated from the indices
    lim_l(l) = (lim_lInd - 1) * bixelWidth + minX - 1 / 2 * bixelWidth;
    lim_r(l) = (lim_rInd - 1) * bixelWidth + minX + 1 / 2 * bixelWidth;
end

% find upmost and downmost leaf pair
topLeafPair = centralLeafPair - maxZ / bixelWidth;
bottomLeafPair = centralLeafPair - minZ / bixelWidth;

% create bool map of active leaf pairs
isActiveLeafPair = zeros(numOfMLCLeafPairs, 1);
isActiveLeafPair(topLeafPair:bottomLeafPair) = 1;

% getting the dimensions of the MLC in order to be able to plot the
% shapes using physical coordinates
MLCWindow = [minX - bixelWidth / 2 maxX + bixelWidth / 2 ...
             minZ - bixelWidth / 2 maxZ + bixelWidth / 2];

geometry.numOfActiveLeafPairs = dimZ;
geometry.leafPairPos          = unique(Z);
geometry.isActiveLeafPair     = isActiveLeafPair;
geometry.centralLeafPair      = centralLeafPair;
geometry.lim_l                = lim_l;
geometry.lim_r                = lim_r;
geometry.bixelIndMap          = bixelIndMap;
geometry.posOfCornerBixel     = [minX 0 minZ];
geometry.MLCWindow            = MLCWindow;

end
