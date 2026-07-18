function updatedInfo = matRad_daoVec2ApertureInfo(apertureInfo, apertureInfoVect)
% matRad function to translate vector representation into struct
% The vector representation of the aperture shape and weight are translated
% into an aperture info struct. At the same time, the updated bixel weight
% vector w is computed and a vector listing the correspondence between leaf
% tips and bixel indices for gradient calculation
%
% call:
%   updatedInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVect)
%
% input:
%   apertureInfo:     aperture shape info struct
%   apertureInfoVect: aperture weights and shapes parameterized as vector
%
% output:
%   updatedInfo: updated aperture shape info struct according to apertureInfoVect
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% function to update the apertureInfo struct after the each iteration of the
% optimization

% This is the static (step-and-shoot) DAO version. VMAT plans must go
% through the matRad_OptimizationProblemVMAT override of this same-named
% static method, which handles arc interpolation, gantry times and the
% analytic bixel-aperture Jacobian.
if apertureInfo.runVMAT
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError(['VMAT apertureInfo passed to the static DAO vector conversion! ' ...
                          'Use matRad_OptimizationProblemVMAT.matRad_daoVec2ApertureInfo (or instance dispatch) instead.']);
end

% initializing variables
updatedInfo = apertureInfo;

updatedInfo.apertureVector = apertureInfoVect;

shapeInd = 1;

indVect = NaN * ones(2 * apertureInfo.doseTotalNumOfLeafPairs, 1);
offset = 0;

% helper function to cope with numerical instabilities through rounding
round2 = @(a, b) round(a * 10^b) / 10^b;

w = zeros(apertureInfo.totalNumOfBixels, 1);

%% update the shapeMaps
% here the new colimator positions are used to create new shapeMaps that
% now include decimal values instead of binary

% loop over all beams
for i = 1:numel(updatedInfo.beam)

    % pre compute left and right bixel edges
    edges_l = updatedInfo.beam(i).posOfCornerBixel(1) + ...
        ((1:size(apertureInfo.beam(i).bixelIndMap, 2)) - 1 - 1 / 2) * updatedInfo.bixelWidth;
    edges_r = updatedInfo.beam(i).posOfCornerBixel(1) + ...
        ((1:size(apertureInfo.beam(i).bixelIndMap, 2)) - 1 + 1 / 2) * updatedInfo.bixelWidth;

    % get dimensions of 2d matrices that store shape/bixel information
    n = apertureInfo.beam(i).numOfActiveLeafPairs;

    % loop over all shapes
    for j = 1:updatedInfo.beam(i).numOfShapes

        % update the shape weight
        updatedInfo.beam(i).shape(j).weight = apertureInfoVect(updatedInfo.beam(i).shape(j).weightOffset) ./ updatedInfo.beam(i).shape(j).jacobiScale;

        % extract left and right leaf positions from shape vector
        vectorIx_L = updatedInfo.beam(i).shape(j).vectorOffset + ((1:n) - 1);
        vectorIx_R = vectorIx_L + apertureInfo.totalNumOfLeafPairs;
        leftLeafPos  = apertureInfoVect(vectorIx_L);
        rightLeafPos = apertureInfoVect(vectorIx_R);

        % update information in shape structure
        updatedInfo.beam(i).shape(j).leftLeafPos  = leftLeafPos;
        updatedInfo.beam(i).shape(j).rightLeafPos = rightLeafPos;

        % rounding for numerical stability
        leftLeafPos  = round2(leftLeafPos, 10);
        rightLeafPos = round2(rightLeafPos, 10);

        % check overshoot of leaf positions
        leftLeafPos(leftLeafPos <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(leftLeafPos <= apertureInfo.beam(i).lim_l);
        rightLeafPos(rightLeafPos <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(rightLeafPos <= apertureInfo.beam(i).lim_l);
        leftLeafPos(leftLeafPos >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(leftLeafPos >= apertureInfo.beam(i).lim_r);
        rightLeafPos(rightLeafPos >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(rightLeafPos >= apertureInfo.beam(i).lim_r);

        %
        xPosIndLeftLeaf  = round((leftLeafPos - apertureInfo.beam(i).posOfCornerBixel(1)) / apertureInfo.bixelWidth + 1);
        xPosIndRightLeaf = round((rightLeafPos - apertureInfo.beam(i).posOfCornerBixel(1)) / apertureInfo.bixelWidth + 1);

        %
        xPosIndLeftLeaf_lim  = floor((apertureInfo.beam(i).lim_l - apertureInfo.beam(i).posOfCornerBixel(1)) / apertureInfo.bixelWidth + 1);
        xPosIndRightLeaf_lim = ceil((apertureInfo.beam(i).lim_r - apertureInfo.beam(i).posOfCornerBixel(1)) / apertureInfo.bixelWidth + 1);

        xPosIndLeftLeaf(xPosIndLeftLeaf <= xPosIndLeftLeaf_lim) = xPosIndLeftLeaf_lim(xPosIndLeftLeaf <= xPosIndLeftLeaf_lim) + 1;
        xPosIndRightLeaf(xPosIndRightLeaf >= xPosIndRightLeaf_lim) = xPosIndRightLeaf_lim(xPosIndRightLeaf >= xPosIndRightLeaf_lim) - 1;

        % check limits because of rounding off issues at maximum, i.e.,
        % enforce round(X.5) -> X
        % LeafPos can occasionally go slightly beyond lim_r, so changed
        % == check to >=
        cornerX = apertureInfo.beam(i).posOfCornerBixel(1);
        overshootL = leftLeafPos >= apertureInfo.beam(i).lim_r;
        overshootR = rightLeafPos >= apertureInfo.beam(i).lim_r;
        xPosIndLeftLeaf(overshootL)  = round(.5 + (leftLeafPos(overshootL) - cornerX) / apertureInfo.bixelWidth);
        xPosIndRightLeaf(overshootR) = round(.5 + (rightLeafPos(overshootR) - cornerX) / apertureInfo.bixelWidth);

        % find the bixel index that the leaves currently touch
        bixelIndMapBeam   = apertureInfo.beam(i).bixelIndMap;
        bixelIndLeftLeaf  = bixelIndMapBeam((xPosIndLeftLeaf - 1) * n + (1:n)');
        bixelIndRightLeaf = bixelIndMapBeam((xPosIndRightLeaf - 1) * n + (1:n)');

        % In non-rectangular fields a leaf tip can land on a column that
        % holds no open bixel (ray) in its row. Snap the touched bixel to
        % the nearest valid (non-NaN) bixel within the row, searching
        % inward (left leaf to the right, right leaf to the left).
        for leafRow = find(isnan(bixelIndLeftLeaf))'
            validCols = find(~isnan(bixelIndMapBeam(leafRow, :)));
            nextCol = validCols(validCols >= xPosIndLeftLeaf(leafRow));
            if ~isempty(nextCol)
                bixelIndLeftLeaf(leafRow) = bixelIndMapBeam(leafRow, nextCol(1));
            end
        end
        for leafRow = find(isnan(bixelIndRightLeaf))'
            validCols = find(~isnan(bixelIndMapBeam(leafRow, :)));
            prevCol = validCols(validCols <= xPosIndRightLeaf(leafRow));
            if ~isempty(prevCol)
                bixelIndRightLeaf(leafRow) = bixelIndMapBeam(leafRow, prevCol(end));
            end
        end

        if any(isnan(bixelIndLeftLeaf)) || any(isnan(bixelIndRightLeaf))
            error('cannot map leaf position to bixel index');
        end

        % store information in index vector for gradient calculation
        indVect(offset + (1:n)) = bixelIndLeftLeaf;
        indVect(offset + (1:n) + apertureInfo.doseTotalNumOfLeafPairs) = bixelIndRightLeaf;
        offset = offset + n;

        % calculate opening fraction for every bixel in shape to construct
        % bixel weight vector

        coveredByLeftLeaf  = bsxfun(@minus, leftLeafPos, edges_l)  / updatedInfo.bixelWidth;
        coveredByRightLeaf = bsxfun(@minus, edges_r, rightLeafPos) / updatedInfo.bixelWidth;

        tempMap = 1 - (coveredByLeftLeaf  + abs(coveredByLeftLeaf))  / 2 - ...
            (coveredByRightLeaf + abs(coveredByRightLeaf)) / 2;

        % find open bixels (ignore positions without a ray in
        % non-rectangular fields, which carry no bixel/dose)
        tempMapIx = tempMap > 0 & ~isnan(apertureInfo.beam(i).bixelIndMap);

        currBixelIx = apertureInfo.beam(i).bixelIndMap(tempMapIx);
        w(currBixelIx) = w(currBixelIx) + tempMap(tempMapIx) * updatedInfo.beam(i).shape(j).weight;

        % save the tempMap (we need to apply a positivity operator !)
        updatedInfo.beam(i).shape(j).shapeMap = (tempMap  + abs(tempMap))  / 2;

        % increment shape index
        shapeInd = shapeInd + 1;

    end
end

% save bixelWeight, apertureVector, and indVect
updatedInfo.bixelWeights = w;
updatedInfo.apertureVector = apertureInfoVect;
updatedInfo.bixelIndices = indVect;

end
