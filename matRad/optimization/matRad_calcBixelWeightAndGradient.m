function accum = matRad_calcBixelWeightAndGradient(calcOptions, mlcOptions, variables, vectorIndices, accum)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate the bixel weights from the aperture vector,
% and also the Jacobian matrix relating these two.
%
% call
%   accum = matRad_calcBixelWeightAndGradient(calcOptions,mlcOptions,variables,vectorIndices,accum)
%
% input
%   calcOptions:    what to compute for this shape (isDAOBeam, saveJacobian)
%   mlcOptions:     MLC geometry of the beam (leaf limits, bixel edges, maps)
%   variables:      leaf positions, weights and times of the current shape
%   vectorIndices:  where this shape's variables live in the aperture vector
%   accum:          accumulators carried across the per-shape calls:
%                     .w             bixel weight vector
%                     .bixelJApVec   Jacobian triplets, fields .vec .i .j
%                     .sumGradSq     squared gradient sum (Jacobi precond.)
%                     .shapeMapW     weighted shape map of the current shape
%                     .counters      running offsets, see .bixelJApVecOffset
%
% output
%   accum:          the same accumulator struct, updated with this shape's
%                   contribution
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

round2 = @(a, b) round(a * 10^b) / 10^b;

%% extract variables from inputs
limLeft = mlcOptions.limLeft;
limRight = mlcOptions.limRight;
edgesLeft = mlcOptions.edgesLeft;
edgesRight = mlcOptions.edgesRight;
centres = mlcOptions.centres;
widths = mlcOptions.widths;
n = mlcOptions.n;
numBix = mlcOptions.numBix;
bixelWidth = mlcOptions.bixelWidth;
bixelIndMap = mlcOptions.bixelIndMap;

weight               = variables.weight;
leftLeafPosInitial   = variables.leftLeafPosInitial;
leftLeafPosFinal     = variables.leftLeafPosFinal;
rightLeafPosInitial  = variables.rightLeafPosInitial;
rightLeafPosFinal    = variables.rightLeafPosFinal;

% Accumulators are unpacked into plain locals here and repacked at the end,
% so that the inner loops below index arrays rather than struct fields.
w              = accum.w;
bixelJApVecVec = accum.bixelJApVec.vec;
bixelJApVecI   = accum.bixelJApVec.i;
bixelJApVecJ   = accum.bixelJApVec.j;
sumGradSq      = accum.sumGradSq;
shapeMapW      = accum.shapeMapW;
counters       = accum.counters;

bixelJApVecOffset = counters.bixelJApVecOffset;

%% sort out order, set up calculation of bixel weight and gradients

% set the initial leaf positions to the minimum leaf positions
% always, instead of the leaf positions at the actual beginning
% of the arc
% this simplifies the calculation
% remember which one is actually I and F in leftMinInd
[leftLeafPosI, leftMinInd] = min([leftLeafPosInitial, leftLeafPosFinal], [], 2);
leftLeafPosF = max([leftLeafPosInitial, leftLeafPosFinal], [], 2);
[rightLeafPosI, rightMinInd] = min([rightLeafPosInitial, rightLeafPosFinal], [], 2);
rightLeafPosF = max([rightLeafPosInitial, rightLeafPosFinal], [], 2);

if calcOptions.saveJacobian
    % only need these variables for the Jacobian

    if calcOptions.isDAOBeam
        jacobiScale = variables.jacobiScale;

        vectorIxLI = vectorIndices.vectorIxLI;
        vectorIxLF = vectorIndices.vectorIxLF;
        vectorIxRI = vectorIndices.vectorIxRI;
        vectorIxRF = vectorIndices.vectorIxRF;
        DAOBeamNumber = vectorIndices.DAOBeamNumber;

        % change the vectorIx_xy elements to remember which
        % apertureVector elements the "new" I and F
        % if leftMinInd is 2, the I and F are switched
        tempL = vectorIxLI;
        tempR = vectorIxRI;
        vectorIxLI(leftMinInd == 2) = vectorIxLF(leftMinInd == 2);
        vectorIxLF(leftMinInd == 2) = tempL(leftMinInd == 2);
        vectorIxRI(rightMinInd == 2) = vectorIxRF(rightMinInd == 2);
        vectorIxRF(rightMinInd == 2) = tempR(rightMinInd == 2);
    else

        weightLast = variables.weightLast;
        weightNext = variables.weightNext;
        jacobiScaleLast = variables.jacobiScaleLast;
        jacobiScaleNext = variables.jacobiScaleNext;

        timeLast = variables.timeLast;
        timeNext = variables.timeNext;
        time = variables.time;

        fracFromLastOptI = variables.fracFromLastOptI;
        fracFromLastOptF = variables.fracFromLastOptF;
        fracFromNextOptI = variables.fracFromNextOptI;
        fracFromNextOptF = variables.fracFromNextOptF;
        fracFromLastOpt = variables.fracFromLastOpt;

        % replicate
        fracFromLastOptI = repmat(fracFromLastOptI, 1, numBix);
        fracFromLastOptF = repmat(fracFromLastOptF, 1, numBix);
        fracFromNextOptI = repmat(fracFromNextOptI, 1, numBix);
        fracFromNextOptF = repmat(fracFromNextOptF, 1, numBix);

        doseAngleBordersDiff = variables.doseAngleBordersDiff;
        doseAngleBordersDiffLast = variables.doseAngleBordersDiffLast;
        doseAngleBordersDiffNext = variables.doseAngleBordersDiffNext;
        timeFactorCurrentLast = variables.timeFactorCurrentLast;
        timeFactorCurrentNext = variables.timeFactorCurrentNext;
        weightFracFromLastDAO = variables.weightFracFromLastDAO;
        timeFracFromLastDAO = variables.timeFracFromLastDAO;
        timeFracFromNextDAO = variables.timeFracFromNextDAO;

        vectorIxLFLast = vectorIndices.vectorIxLFLast;
        vectorIxLINext = vectorIndices.vectorIxLINext;
        vectorIxRFLast = vectorIndices.vectorIxRFLast;
        vectorIxRINext = vectorIndices.vectorIxRINext;
        DAOBeamNumberLast   = vectorIndices.DAOBeamNumberLast;
        DAOBeamNumberNext   = vectorIndices.DAOBeamNumberNext;
        tIxLast        = vectorIndices.tIxLast;
        tIxNext        = vectorIndices.tIxNext;

        tempL = vectorIxLFLast;
        tempR = vectorIxRFLast;

        vectorIxLFLast(leftMinInd == 2) = vectorIxLINext(leftMinInd == 2);
        vectorIxLINext(leftMinInd == 2) = tempL(leftMinInd == 2);

        vectorIxRFLast(rightMinInd == 2) = vectorIxRINext(rightMinInd == 2);
        vectorIxRINext(rightMinInd == 2) = tempR(rightMinInd == 2);
    end
end

leftLeafPosI = round2(leftLeafPosI, 10);
leftLeafPosF = round2(leftLeafPosF, 10);
rightLeafPosI = round2(rightLeafPosI, 10);
rightLeafPosF = round2(rightLeafPosF, 10);

leftLeafPosI(leftLeafPosI <= limLeft) = limLeft(leftLeafPosI <= limLeft);
leftLeafPosF(leftLeafPosF <= limLeft) = limLeft(leftLeafPosF <= limLeft);
rightLeafPosI(rightLeafPosI <= limLeft) = limLeft(rightLeafPosI <= limLeft);
rightLeafPosF(rightLeafPosF <= limLeft) = limLeft(rightLeafPosF <= limLeft);
leftLeafPosI(leftLeafPosI >= limRight) = limRight(leftLeafPosI >= limRight);
leftLeafPosF(leftLeafPosF >= limRight) = limRight(leftLeafPosF >= limRight);
rightLeafPosI(rightLeafPosI >= limRight) = limRight(rightLeafPosI >= limRight);
rightLeafPosF(rightLeafPosF >= limRight) = limRight(rightLeafPosF >= limRight);

% find bixel indices where leaves are located
xPosIndLeftLeafI = min(floor((leftLeafPosI - edgesLeft(1)) ./ bixelWidth) + 1, numBix);
xPosIndLeftLeafF = min(floor((leftLeafPosF - edgesLeft(1)) ./ bixelWidth) + 1, numBix);
xPosIndRightLeafI = min(floor((rightLeafPosI - edgesLeft(1)) ./ bixelWidth) + 1, numBix);
xPosIndRightLeafF = min(floor((rightLeafPosF - edgesLeft(1)) ./ bixelWidth) + 1, numBix);
%
xPosLinearIndLeftLeafI = sub2ind([n numBix], (1:n)', xPosIndLeftLeafI);
xPosLinearIndLeftLeafF = sub2ind([n numBix], (1:n)', xPosIndLeftLeafF);
xPosLinearIndRightLeafI = sub2ind([n numBix], (1:n)', xPosIndRightLeafI);
xPosLinearIndRightLeafF = sub2ind([n numBix], (1:n)', xPosIndRightLeafF);

% distance each leaf sweeps, and the bixel edges/widths at the bixels the
% initial and final leaf positions fall into. Hoisted out of the overshoot
% corrections and gradients below, which all reuse them.
leftLeafSpan = leftLeafPosF - leftLeafPosI;
rightLeafSpan = rightLeafPosF - rightLeafPosI;
edgeLeftAtLI = edgesLeft(xPosIndLeftLeafI)';
edgeRightAtLF = edgesRight(xPosIndLeftLeafF)';
edgeLeftAtRI = edgesLeft(xPosIndRightLeafI)';
edgeRightAtRF = edgesRight(xPosIndRightLeafF)';
widthAtLI = widths(xPosIndLeftLeafI)';
widthAtLF = widths(xPosIndLeftLeafF)';
widthAtRI = widths(xPosIndRightLeafI)';
widthAtRF = widths(xPosIndRightLeafF)';

%
% leaves sweep from _I to _F, with weight
%

%% bixel weight calculation

% calculate fraction of fluence uncovered by left leaf
% initial computation
uncoveredByLeftLeaf = bsxfun(@minus, centres, leftLeafPosI) ./ repmat(leftLeafSpan, 1, numBix);
% correct for overshoot in initial and final leaf positions
uncoveredByLeftLeaf(xPosLinearIndLeftLeafI) = uncoveredByLeftLeaf(xPosLinearIndLeftLeafI) + ...
    (leftLeafPosI - edgeLeftAtLI).^2 ./ (leftLeafSpan .* widthAtLI .* 2);
uncoveredByLeftLeaf(xPosLinearIndLeftLeafF) = uncoveredByLeftLeaf(xPosLinearIndLeftLeafF) - ...
    (edgeRightAtLF - leftLeafPosF).^2 ./ (leftLeafSpan .* widthAtLF .* 2);
% round <0 to 0, >1 to 1
uncoveredByLeftLeaf(uncoveredByLeftLeaf < 0) = 0;
uncoveredByLeftLeaf(uncoveredByLeftLeaf > 1) = 1;

% calculate fraction of fluence covered by right leaf
% initial computation
coveredByRightLeaf = bsxfun(@minus, centres, rightLeafPosI) ./ repmat(rightLeafSpan, 1, numBix);
% correct for overshoot in initial and final leaf positions
coveredByRightLeaf(xPosLinearIndRightLeafI) = coveredByRightLeaf(xPosLinearIndRightLeafI) + ...
    (rightLeafPosI - edgeLeftAtRI).^2 ./ (rightLeafSpan .* widthAtRI .* 2);
coveredByRightLeaf(xPosLinearIndRightLeafF) = coveredByRightLeaf(xPosLinearIndRightLeafF) - ...
    (edgeRightAtRF - rightLeafPosF).^2 ./ (rightLeafSpan .* widthAtRF .* 2);
% round <0 to 0, >1 to 1
coveredByRightLeaf(coveredByRightLeaf < 0) = 0;
coveredByRightLeaf(coveredByRightLeaf > 1) = 1;

%% gradient calculation

dUl_dLI = bsxfun(@minus, centres, leftLeafPosF) ./ (repmat(leftLeafSpan, 1, numBix)).^2;
dUl_dLF = bsxfun(@minus, leftLeafPosI, centres) ./ (repmat(leftLeafSpan, 1, numBix)).^2;

dCr_dRI = bsxfun(@minus, centres, rightLeafPosF) ./ (repmat(rightLeafSpan, 1, numBix)).^2;
dCr_dRF = bsxfun(@minus, rightLeafPosI, centres) ./ (repmat(rightLeafSpan, 1, numBix)).^2;

dUl_dLI(xPosLinearIndLeftLeafI) = dUl_dLI(xPosLinearIndLeftLeafI) + ...
    ((leftLeafPosI - edgeLeftAtLI) .* (2 * leftLeafPosF - leftLeafPosI - edgeLeftAtLI)) ./ (leftLeafSpan.^2 .* widthAtLI .* 2);
dUl_dLF(xPosLinearIndLeftLeafI) = dUl_dLF(xPosLinearIndLeftLeafI) - ((leftLeafPosI - edgeLeftAtLI).^2) ./ (leftLeafSpan.^2 .* widthAtLI .* 2);
dUl_dLI(xPosLinearIndLeftLeafF) = dUl_dLI(xPosLinearIndLeftLeafF) - ((edgeRightAtLF - leftLeafPosF).^2) ./ (leftLeafSpan.^2 .* widthAtLF .* 2);
dUl_dLF(xPosLinearIndLeftLeafF) = dUl_dLF(xPosLinearIndLeftLeafF) + ...
    ((edgeRightAtLF - leftLeafPosF) .* (leftLeafPosF + edgeRightAtLF - 2 * leftLeafPosI)) ./ (leftLeafSpan.^2 .* widthAtLF .* 2);

dCr_dRI(xPosLinearIndRightLeafI) = dCr_dRI(xPosLinearIndRightLeafI) + ...
    ((rightLeafPosI - edgeLeftAtRI) .* (2 * rightLeafPosF - rightLeafPosI - edgeLeftAtRI)) ./ (rightLeafSpan.^2 .* widthAtRI .* 2);
dCr_dRF(xPosLinearIndRightLeafI) = dCr_dRF(xPosLinearIndRightLeafI) - ((rightLeafPosI - edgeLeftAtRI).^2) ./ (rightLeafSpan.^2 .* widthAtRI .* 2);
dCr_dRI(xPosLinearIndRightLeafF) = dCr_dRI(xPosLinearIndRightLeafF) - ((edgeRightAtRF - rightLeafPosF).^2) ./ (rightLeafSpan.^2 .* widthAtRF .* 2);
dCr_dRF(xPosLinearIndRightLeafF) = dCr_dRF(xPosLinearIndRightLeafF) + ...
    ((edgeRightAtRF - rightLeafPosF) .* (rightLeafPosF + edgeRightAtRF - 2 * rightLeafPosI)) ./ (rightLeafSpan.^2 .* widthAtRF .* 2);

for k = 1:n
    dUl_dLI(k, 1:(xPosIndLeftLeafI(k) - 1)) = 0;
    dUl_dLF(k, 1:(xPosIndLeftLeafI(k) - 1)) = 0;
    dUl_dLI(k, (xPosIndLeftLeafF(k) + 1):numBix) = 0;
    dUl_dLF(k, (xPosIndLeftLeafF(k) + 1):numBix) = 0;

    if xPosIndLeftLeafI(k) >= xPosIndLeftLeafF(k)
        % in discrete aperture, the xPosIndLeftLeafI is greater than
        % xPosIndLeftLeafM when leaf positions are at a bixel boundary

        % 19 July 2017 in journal
        dUl_dLI(k, xPosIndLeftLeafI(k)) = -1 / (2 * widths(xPosIndLeftLeafI(k))');
        dUl_dLF(k, xPosIndLeftLeafF(k)) = -1 / (2 * widths(xPosIndLeftLeafF(k))');
        if leftLeafPosF(k) - leftLeafPosI(k) <= eps(max(limRight))
            uncoveredByLeftLeaf(k, xPosIndLeftLeafI(k)) = (edgesRight(xPosIndLeftLeafI(k)) - leftLeafPosI(k)) ./ widths(xPosIndLeftLeafI(k));
            uncoveredByLeftLeaf(k, xPosIndLeftLeafF(k)) = (edgesRight(xPosIndLeftLeafF(k)) - leftLeafPosF(k)) ./ widths(xPosIndLeftLeafF(k));
        end
    end

    dCr_dRI(k, 1:(xPosIndRightLeafI(k) - 1)) = 0;
    dCr_dRF(k, 1:(xPosIndRightLeafI(k) - 1)) = 0;
    dCr_dRI(k, (xPosIndRightLeafF(k) + 1):numBix) = 0;
    dCr_dRF(k, (xPosIndRightLeafF(k) + 1):numBix) = 0;

    if xPosIndRightLeafI(k) >= xPosIndRightLeafF(k)
        dCr_dRI(k, xPosIndRightLeafI(k)) = -1 / (2 * widths(xPosIndRightLeafI(k))');
        dCr_dRF(k, xPosIndRightLeafF(k)) = -1 / (2 * widths(xPosIndRightLeafF(k))');
        if rightLeafPosF(k) - rightLeafPosI(k) <= eps(max(limRight))
            coveredByRightLeaf(k, xPosIndRightLeafI(k)) = (edgesRight(xPosIndRightLeafI(k)) - rightLeafPosI(k)) ./ widths(xPosIndRightLeafI(k));
            coveredByRightLeaf(k, xPosIndRightLeafF(k)) = (edgesRight(xPosIndRightLeafF(k)) - rightLeafPosF(k)) ./ widths(xPosIndRightLeafF(k));
        end
    end
end

% store information for Jacobi preconditioning
sumGradSq = sumGradSq + ...
    mean([sum((dUl_dLI).^2, 2); sum((dUl_dLF).^2, 2); sum((dUl_dLF).^2, 2); sum((dCr_dRI).^2, 2); sum((dCr_dRF).^2, 2); sum((dCr_dRF).^2, 2)]);

%% save the bixel weights
% fluence is equal to fluence not covered by left leaf minus
% fluence covered by left leaf
shapeMap = uncoveredByLeftLeaf - coveredByRightLeaf;
shapeMap = round2(shapeMap, 15);
shapeMap(isnan(shapeMap)) = 0;

% find open bixels
% shapeMapIx = shapeMap > 0;
shapeMapIx = ~isnan(bixelIndMap);

currBixelIx = bixelIndMap(shapeMapIx);
w(currBixelIx) = w(currBixelIx) + shapeMap(shapeMapIx) .* weight;
shapeMapW = shapeMapW + shapeMap .* weight;

%% save the gradients

if calcOptions.saveJacobian

    numSaveBixel = nnz(shapeMapIx);

    if calcOptions.isDAOBeam
        % indices
        vectorIxMatLI = repmat(vectorIxLI', 1, numBix);
        vectorIxMatLF = repmat(vectorIxLF', 1, numBix);
        vectorIxMatRI = repmat(vectorIxRI', 1, numBix);
        vectorIxMatRF = repmat(vectorIxRF', 1, numBix);

        % wrt weight
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = shapeMap(shapeMapIx) ./ jacobiScale;
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = DAOBeamNumber;
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;

        % wrt initial left
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = dUl_dLI(shapeMapIx) .* weight;
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = vectorIxMatLI(shapeMapIx);
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;

        % wrt final left
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = dUl_dLF(shapeMapIx) .* weight;
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = vectorIxMatLF(shapeMapIx);
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;

        % wrt initial right
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = -dCr_dRI(shapeMapIx) .* weight;
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = vectorIxMatRI(shapeMapIx);
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;

        % wrt final right
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = -dCr_dRF(shapeMapIx) .* weight;
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = vectorIxMatRF(shapeMapIx);
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;

    else
        % indices
        vectorIxMatLFLast = repmat(vectorIxLFLast', 1, numBix);
        vectorIxMatLINext = repmat(vectorIxLINext', 1, numBix);
        vectorIxMatRFLast = repmat(vectorIxRFLast', 1, numBix);
        vectorIxMatRINext = repmat(vectorIxRINext', 1, numBix);

        % wrt last weight
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = fracFromLastOpt * (time ./ timeLast) * shapeMap(shapeMapIx) ./ jacobiScaleLast;
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = DAOBeamNumberLast;
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;

        % wrt next weight
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = (1 - fracFromLastOpt) * (time ./ timeNext) * shapeMap(shapeMapIx) ./ jacobiScaleNext;
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = DAOBeamNumberNext;
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;

        % wrt initial left (optimization vector)
        % initial (interpolated arc)
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = fracFromLastOptI(shapeMapIx) .* dUl_dLI(shapeMapIx) .* weight;
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = vectorIxMatLFLast(shapeMapIx);
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;
        % final (interpolated arc)
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = fracFromLastOptF(shapeMapIx) .* dUl_dLF(shapeMapIx) .* weight;
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = vectorIxMatLFLast(shapeMapIx);
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;

        % wrt final left (optimization vector)
        % initial (interpolated arc)
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = fracFromNextOptI(shapeMapIx) .* dUl_dLI(shapeMapIx) .* weight;
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = vectorIxMatLINext(shapeMapIx);
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;
        % final (interpolated arc)
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = fracFromNextOptF(shapeMapIx) .* dUl_dLF(shapeMapIx) .* weight;
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = vectorIxMatLINext(shapeMapIx);
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;

        % wrt initial right (optimization vector)
        % initial (interpolated arc)
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = -fracFromLastOptI(shapeMapIx) .* dCr_dRI(shapeMapIx) .* weight;
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = vectorIxMatRFLast(shapeMapIx);
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;
        % final (interpolated arc)
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = -fracFromLastOptF(shapeMapIx) .* dCr_dRF(shapeMapIx) .* weight;
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = vectorIxMatRFLast(shapeMapIx);
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;

        % wrt final right (optimization vector)
        % initial (interpolated arc)
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = -fracFromNextOptI(shapeMapIx) .* dCr_dRI(shapeMapIx) .* weight;
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = vectorIxMatRINext(shapeMapIx);
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;
        % final (interpolated arc)
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = -fracFromNextOptF(shapeMapIx) .* dCr_dRF(shapeMapIx) .* weight;
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = vectorIxMatRINext(shapeMapIx);
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;

        % wrt last time
        dTimeLast = -weightFracFromLastDAO .* timeFracFromNextDAO .* ...
                    (weightLast ./ doseAngleBordersDiffNext) .* (timeNext ./ timeLast.^2) + ...
                    (1 - weightFracFromLastDAO) .* timeFracFromLastDAO .* ...
                    (weightNext ./ doseAngleBordersDiffLast) .* (1 ./ timeNext);
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = ...
            (doseAngleBordersDiff .* timeFactorCurrentLast) .* dTimeLast * shapeMap(shapeMapIx);
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = tIxLast;
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;

        % wrt next time
        dTimeNext = weightFracFromLastDAO .* timeFracFromNextDAO .* ...
                    (weightLast ./ doseAngleBordersDiffNext) .* (1 ./ timeLast) - ...
                    (1 - weightFracFromLastDAO) .* timeFracFromLastDAO .* ...
                    (weightNext ./ doseAngleBordersDiffLast) .* (timeLast ./ timeNext.^2);
        bixelJApVecVec(bixelJApVecOffset + (1:numSaveBixel)) = ...
            (doseAngleBordersDiff .* timeFactorCurrentNext) .* dTimeNext * shapeMap(shapeMapIx);
        bixelJApVecI(bixelJApVecOffset + (1:numSaveBixel)) = tIxNext;
        bixelJApVecJ(bixelJApVecOffset + (1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVecOffset = bixelJApVecOffset + numSaveBixel;
    end

end

% update counters and repack the accumulators
counters.bixelJApVecOffset = bixelJApVecOffset;

accum.w               = w;
accum.bixelJApVec.vec = bixelJApVecVec;
accum.bixelJApVec.i   = bixelJApVecI;
accum.bixelJApVec.j   = bixelJApVecJ;
accum.sumGradSq       = sumGradSq;
accum.shapeMapW       = shapeMapW;
accum.counters        = counters;

end
