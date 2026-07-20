function updatedInfo = matRad_recalcApertureBixelWeights(apertureInfo)
% matRad function to compute bixel weights for a recalculated (fine-angle)
% VMAT aperture sequence. Each beam's dose sector is integrated as two
% uniform leaf sweeps (initial position -> centre, centre -> final
% position) with the corresponding half-sector weights, reusing the
% dynamic-sweep fluence model of matRad_bixWeightAndGrad (no Jacobian).
%
% call
%   updatedInfo = matRad_recalcApertureBixelWeights(apertureInfo)
%
% input
%   apertureInfo:  recalculated aperture info struct (all beams carry
%                  shape(1) with weight/weight_I/weight_F and leaf
%                  positions; see matRad_recalcApertureInfo)
%
% output
%   updatedInfo:   apertureInfo with updated shapeMaps and bixelWeights
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

updatedInfo = apertureInfo;

% options for the bixel weight calculation - forward fluence only
mlcOptions.bixelWidth = apertureInfo.bixelWidth;
calcOptions.continuousAperture = updatedInfo.continuousAperture;
calcOptions.saveJacobian = false;
vectorIndices.totalNumOfShapes = apertureInfo.totalNumOfShapes;

w = zeros(apertureInfo.totalNumOfBixels, 1);

% dummy jacobian containers (unused with saveJacobian = false)
bixelJApVecVec = [];
bixelJApVecI = 0;
bixelJApVecJ = 0;
counters.bixelJApVecOffset = 0;

% loop over all beams
for i = 1:numel(updatedInfo.beam)

    % pre compute left and right bixel edges
    edgesLeft = updatedInfo.beam(i).posOfCornerBixel(1) + ...
        ((1:size(apertureInfo.beam(i).bixelIndMap, 2)) - 1 - 1 / 2) * updatedInfo.bixelWidth;
    edgesRight = updatedInfo.beam(i).posOfCornerBixel(1) + ...
        ((1:size(apertureInfo.beam(i).bixelIndMap, 2)) - 1 + 1 / 2) * updatedInfo.bixelWidth;

    n = apertureInfo.beam(i).numOfActiveLeafPairs;

    mlcOptions.limLeft = apertureInfo.beam(i).limLeft;
    mlcOptions.limRight = apertureInfo.beam(i).limRight;
    mlcOptions.edgesLeft = edgesLeft;
    mlcOptions.edgesRight = edgesRight;
    mlcOptions.centres = (edgesLeft + edgesRight) / 2;
    mlcOptions.widths = edgesRight - edgesLeft;
    mlcOptions.n = n;
    mlcOptions.numBix = size(apertureInfo.beam(i).bixelIndMap, 2);
    mlcOptions.bixelIndMap = apertureInfo.beam(i).bixelIndMap;
    calcOptions.isDAOBeam = updatedInfo.arc.beam(i).isDAOBeam;

    shapeMap_I = zeros(size(updatedInfo.beam(i).bixelIndMap));
    shapeMap_F = zeros(size(updatedInfo.beam(i).bixelIndMap));
    sumGradSq = 0;

    % extract the weights and leaf positions from the apertureInfo
    weight = updatedInfo.beam(i).shape(1).weight;
    if isfield(updatedInfo.beam(i).shape(1), 'weight_I')
        weight_I = updatedInfo.beam(i).shape(1).weight_I;
        weight_F = updatedInfo.beam(i).shape(1).weight_F;
    else
        % only happens at original angular resolution
        weight_I = weight .* updatedInfo.arc.beam(i).doseAngleBorderCentreDiff(1) ./ updatedInfo.arc.beam(i).doseAngleBordersDiff;
        weight_F = weight .* updatedInfo.arc.beam(i).doseAngleBorderCentreDiff(2) ./ updatedInfo.arc.beam(i).doseAngleBordersDiff;
    end

    if weight_I + weight_F ~= weight
        % sometimes the sum is different than one by ~10^-16
        % (rounding error in the division)
        weight_F = weight - weight_I;
    end

    % initial half sector: sweep from the initial to the central leaf
    % position with the initial half-sector weight
    variables.weight = weight_I;
    if updatedInfo.continuousAperture
        variables.leftLeafPosInitial     = updatedInfo.beam(i).shape(1).leftLeafPosInitial;
        variables.leftLeafPosFinal     = updatedInfo.beam(i).shape(1).leftLeafPos;
        variables.rightLeafPosInitial    = updatedInfo.beam(i).shape(1).rightLeafPosInitial;
        variables.rightLeafPosFinal    = updatedInfo.beam(i).shape(1).rightLeafPos;
    else
        variables.leftLeafPosInitial     = updatedInfo.beam(i).shape(1).leftLeafPos;
        variables.leftLeafPosFinal     = updatedInfo.beam(i).shape(1).leftLeafPos;
        variables.rightLeafPosInitial    = updatedInfo.beam(i).shape(1).rightLeafPos;
        variables.rightLeafPosFinal    = updatedInfo.beam(i).shape(1).rightLeafPos;
    end

    accum.w               = w;
    accum.bixelJApVec.vec = bixelJApVecVec;
    accum.bixelJApVec.i   = bixelJApVecI;
    accum.bixelJApVec.j   = bixelJApVecJ;
    accum.sumGradSq       = sumGradSq;
    accum.shapeMapW       = shapeMap_I;
    accum.counters        = counters;

    accum = matRad_bixWeightAndGrad(calcOptions, mlcOptions, variables, vectorIndices, accum);

    w          = accum.w;
    sumGradSq  = accum.sumGradSq;
    shapeMap_I = accum.shapeMapW;
    counters   = accum.counters;

    % final half sector: sweep from the central to the final leaf
    % position with the final half-sector weight
    variables.weight = weight_F;
    if updatedInfo.continuousAperture
        variables.leftLeafPosInitial     = updatedInfo.beam(i).shape(1).leftLeafPos;
        variables.leftLeafPosFinal     = updatedInfo.beam(i).shape(1).leftLeafPosFinal;
        variables.rightLeafPosInitial    = updatedInfo.beam(i).shape(1).rightLeafPos;
        variables.rightLeafPosFinal    = updatedInfo.beam(i).shape(1).rightLeafPosFinal;
    end

    accum.w               = w;
    accum.bixelJApVec.vec = bixelJApVecVec;
    accum.bixelJApVec.i   = bixelJApVecI;
    accum.bixelJApVec.j   = bixelJApVecJ;
    accum.sumGradSq       = sumGradSq;
    accum.shapeMapW       = shapeMap_F;
    accum.counters        = counters;

    accum = matRad_bixWeightAndGrad(calcOptions, mlcOptions, variables, vectorIndices, accum);

    w          = accum.w;
    shapeMap_F = accum.shapeMapW;
    counters   = accum.counters;

    updatedInfo.beam(i).shape(1).shapeMap = shapeMap_I + shapeMap_F;

end

updatedInfo.bixelWeights = w;

end
