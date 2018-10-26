function [w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j,sumGradSq,shapeMapW,counters] = ...
    matRad_bixWeightAndGrad(calcOptions,mlcOptions,variables,vectorIndices,counters,w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j,sumGradSq,shapeMapW)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate the bixel weights from the aperture vector,
% and also the Jacobian matrix relating these two.
%
% call
%   [w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j,sqrtSumGradSq,shapeMap,counters] = ...
%    matRad_bixWeightAndGrad(calcOptions,mlcOptions,variables,vectorIndices,counters,w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j)
%
% input
%
% output
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018, Mark Bangert, on behalf of the matRad development team
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

round2 = @(a,b) round(a*10^b)/10^b;

%% extract variables from inputs
lim_l = mlcOptions.lim_l;
lim_r = mlcOptions.lim_r;
edges_l = mlcOptions.edges_l;
edges_r = mlcOptions.edges_r;
centres = mlcOptions.centres;
widths = mlcOptions.widths;
n = mlcOptions.n;
numBix = mlcOptions.numBix;
bixelWidth = mlcOptions.bixelWidth;
bixelIndMap = mlcOptions.bixelIndMap;

weight_I        = variables.weight_I;
weight_F        = variables.weight_F;
weightFactor_I  = variables.weightFactor_I;
weightFactor_F  = variables.weightFactor_F;
leftLeafPos_I   = variables.leftLeafPos_I;
leftLeafPos_F   = variables.leftLeafPos_F;
rightLeafPos_I  = variables.rightLeafPos_I;
rightLeafPos_F  = variables.rightLeafPos_F;

bixelJApVec_offset = counters.bixelJApVec_offset;

totalNumOfShapes = vectorIndices.totalNumOfShapes;


%% sort out order, set up calculation of bixel weight and gradients

% set the initial leaf positions to the minimum leaf positions
% always, instead of the leaf positions at the actual beginning
% of the arc
% this simplifies the calculation
% remember which one is actually I and F in leftMinInd
[leftLeafPosI,leftMinInd] = min([leftLeafPos_I,leftLeafPos_F],[],2);
leftLeafPosF = max([leftLeafPos_I,leftLeafPos_F],[],2);
[rightLeafPosI,rightMinInd] = min([rightLeafPos_I,rightLeafPos_F],[],2);
rightLeafPosF = max([rightLeafPos_I,rightLeafPos_F],[],2);

if calcOptions.saveJacobian
    % only need these variables for the Jacobian
    
    probability_dTVec   = variables.probability_dTVec;
    tIx_Vec             = vectorIndices.tIx_Vec;
    
    if calcOptions.DAOBeam
        jacobiScale_I = variables.jacobiScale_I;
        jacobiScale_F = variables.jacobiScale_F;
        
        vectorIx_LI = vectorIndices.vectorIx_LI;
        vectorIx_LF = vectorIndices.vectorIx_LF;
        vectorIx_RI = vectorIndices.vectorIx_RI;
        vectorIx_RF = vectorIndices.vectorIx_RF;
        DAOindex = vectorIndices.DAOindex;
        
        % change the vectorIx_xy elements to remember which
        % apertureVector elements the "new" I and F
        % if leftMinInd is 2, the I and F are switched
        tempL = vectorIx_LI;
        tempR = vectorIx_RI;
        vectorIx_LI(leftMinInd == 2) = vectorIx_LF(leftMinInd == 2);
        vectorIx_LF(leftMinInd == 2) = tempL(leftMinInd == 2);
        vectorIx_RI(rightMinInd == 2) = vectorIx_RF(rightMinInd == 2);
        vectorIx_RF(rightMinInd == 2) = tempR(rightMinInd == 2);
    else
        
        weight_last_I = variables.weight_last_I;
        weight_last_F = variables.weight_last_F;
        weight_next_I = variables.weight_next_I;
        weight_next_F = variables.weight_next_F;
        jacobiScale_last_I = variables.jacobiScale_last_I;
        jacobiScale_last_F = variables.jacobiScale_last_F;
        jacobiScale_next_I = variables.jacobiScale_next_I;
        jacobiScale_next_F = variables.jacobiScale_next_F;
        
        time_last = variables.time_last;
        time_next = variables.time_next;
        time = variables.time;
        
        fracFromLastOptI = variables.fracFromLastOptI;
        fracFromLastOptF = variables.fracFromLastOptF;
        fracFromNextOptI = variables.fracFromNextOptI;
        fracFromNextOptF = variables.fracFromNextOptF;
        fracFromLastOpt = variables.fracFromLastOpt;
        
        doseAngleBordersDiff = variables.doseAngleBordersDiff;
        doseAngleBordersDiff_last = variables.doseAngleBordersDiff_last;
        doseAngleBordersDiff_next = variables.doseAngleBordersDiff_next;
        timeFacCurr_last = variables.timeFacCurr_last;
        timeFacCurr_next = variables.timeFacCurr_next;
        fracFromLastDAO = variables.fracFromLastDAO;
        timeFracFromLastDAO = variables.timeFracFromLastDAO;
        timeFracFromNextDAO = variables.timeFracFromNextDAO;
        
        vectorIx_LF_last = vectorIndices.vectorIx_LF_last;
        vectorIx_LI_next = vectorIndices.vectorIx_LI_next;
        vectorIx_RF_last = vectorIndices.vectorIx_RF_last;
        vectorIx_RI_next = vectorIndices.vectorIx_RI_next;
        DAOindex_last   = vectorIndices.DAOindex_last;
        DAOindex_next   = vectorIndices.DAOindex_next;
        tIx_last        = vectorIndices.tIx_last;
        tIx_next        = vectorIndices.tIx_next;
        
        tempL = vectorIx_LF_last;
        tempR = vectorIx_RF_last;
        
        vectorIx_LF_last(leftMinInd == 2) = vectorIx_LI_next(leftMinInd == 2);
        vectorIx_LI_next(leftMinInd == 2) = tempL(leftMinInd == 2);
        
        vectorIx_RF_last(rightMinInd == 2) = vectorIx_RI_next(rightMinInd == 2);
        vectorIx_RI_next(rightMinInd == 2) = tempR(rightMinInd == 2);
    end
end

leftLeafPosI = round2(leftLeafPosI,10);
leftLeafPosF = round2(leftLeafPosF,10);
rightLeafPosI = round2(rightLeafPosI,10);
rightLeafPosF = round2(rightLeafPosF,10);

leftLeafPosI(leftLeafPosI <= lim_l) = lim_l(leftLeafPosI <= lim_l);
leftLeafPosF(leftLeafPosF <= lim_l) = lim_l(leftLeafPosF <= lim_l);
rightLeafPosI(rightLeafPosI <= lim_l) = lim_l(rightLeafPosI <= lim_l);
rightLeafPosF(rightLeafPosF <= lim_l) = lim_l(rightLeafPosF <= lim_l);
leftLeafPosI(leftLeafPosI >= lim_r) = lim_r(leftLeafPosI >= lim_r);
leftLeafPosF(leftLeafPosF >= lim_r) = lim_r(leftLeafPosF >= lim_r);
rightLeafPosI(rightLeafPosI >= lim_r) = lim_r(rightLeafPosI >= lim_r);
rightLeafPosF(rightLeafPosF >= lim_r) = lim_r(rightLeafPosF >= lim_r);

% determine middle leaf positions
leftLeafPosM   = weightFactor_F.*leftLeafPosI+weightFactor_I.*leftLeafPosF;
rightLeafPosM  = weightFactor_F.*rightLeafPosI+weightFactor_I.*rightLeafPosF;

% find bixel indices where leaves are located
xPosIndLeftLeafI = min(floor((leftLeafPosI-edges_l(1))./bixelWidth)+1,numBix);
xPosIndLeftLeafM = min(floor((leftLeafPosM-edges_l(1))./bixelWidth)+1,numBix);
xPosIndLeftLeafF = max(ceil((leftLeafPosF-edges_r(1))./bixelWidth)+1,1);
xPosIndRightLeafI = min(floor((rightLeafPosI-edges_l(1))./bixelWidth)+1,numBix);
xPosIndRightLeafM = min(floor((rightLeafPosM-edges_l(1))./bixelWidth)+1,numBix);
xPosIndRightLeafF = max(ceil((rightLeafPosF-edges_r(1))./bixelWidth)+1,1);
%
xPosLinearIndLeftLeafI = sub2ind([n numBix],(1:n)',xPosIndLeftLeafI);
xPosLinearIndLeftLeafM = sub2ind([n numBix],(1:n)',xPosIndLeftLeafM);
xPosLinearIndLeftLeafF = sub2ind([n numBix],(1:n)',xPosIndLeftLeafF);
xPosLinearIndRightLeafI = sub2ind([n numBix],(1:n)',xPosIndRightLeafI);
xPosLinearIndRightLeafM = sub2ind([n numBix],(1:n)',xPosIndRightLeafM);
xPosLinearIndRightLeafF = sub2ind([n numBix],(1:n)',xPosIndRightLeafF);


%
% leaves sweep from _I to _M, with weight weight_I.*weightFactor_I
%


%% bixel weight calculation

%calculate fraction of fluence uncovered by left leaf
%initial computation
uncoveredByLeftLeaf = bsxfun(@minus,centres,leftLeafPosI)./repmat(leftLeafPosM-leftLeafPosI,1,numBix);
%correct for overshoot in initial and final leaf positions
uncoveredByLeftLeaf(xPosLinearIndLeftLeafI) = uncoveredByLeftLeaf(xPosLinearIndLeftLeafI) + (leftLeafPosI-edges_l(xPosIndLeftLeafI)').^2./((leftLeafPosM-leftLeafPosI).*(widths(xPosIndLeftLeafI)').*2);
uncoveredByLeftLeaf(xPosLinearIndLeftLeafM) = uncoveredByLeftLeaf(xPosLinearIndLeftLeafM) - (edges_r(xPosIndLeftLeafM)'-leftLeafPosM).^2./((leftLeafPosM-leftLeafPosI).*(widths(xPosIndLeftLeafM)').*2);
%round <0 to 0, >1 to 1
uncoveredByLeftLeaf(uncoveredByLeftLeaf < 0) = 0;
uncoveredByLeftLeaf(uncoveredByLeftLeaf > 1) = 1;

%calculate fraction of fluence covered by right leaf
%initial computation
coveredByRightLeaf = bsxfun(@minus,centres,rightLeafPosI)./repmat(rightLeafPosM-rightLeafPosI,1,numBix);
%correct for overshoot in initial and final leaf positions
coveredByRightLeaf(xPosLinearIndRightLeafI) = coveredByRightLeaf(xPosLinearIndRightLeafI) + (rightLeafPosI-edges_l(xPosIndRightLeafI)').^2./((rightLeafPosM-rightLeafPosI).*(widths(xPosIndRightLeafI)').*2);
coveredByRightLeaf(xPosLinearIndRightLeafM) = coveredByRightLeaf(xPosLinearIndRightLeafM) - (edges_r(xPosIndRightLeafM)'-rightLeafPosM).^2./((rightLeafPosM-rightLeafPosI).*(widths(xPosIndRightLeafM)').*2);
%round <0 to 0, >1 to 1
coveredByRightLeaf(coveredByRightLeaf < 0) = 0;
coveredByRightLeaf(coveredByRightLeaf > 1) = 1;

%% gradient calculation

dUl_dLI = bsxfun(@minus,centres,leftLeafPosM)./(repmat(leftLeafPosM-leftLeafPosI,1,numBix)).^2;
dUl_dLM = bsxfun(@minus,leftLeafPosI,centres)./(repmat(leftLeafPosM-leftLeafPosI,1,numBix)).^2;

dCr_dRI = bsxfun(@minus,centres,rightLeafPosM)./(repmat(rightLeafPosM-rightLeafPosI,1,numBix)).^2;
dCr_dRM = bsxfun(@minus,rightLeafPosI,centres)./(repmat(rightLeafPosM-rightLeafPosI,1,numBix)).^2;

dUl_dLI(xPosLinearIndLeftLeafI) = dUl_dLI(xPosLinearIndLeftLeafI) + ((leftLeafPosI-edges_l(xPosIndLeftLeafI)').*(2*leftLeafPosM-leftLeafPosI-edges_l(xPosIndLeftLeafI)'))./((leftLeafPosM-leftLeafPosI).^2.*(widths(xPosIndLeftLeafI)').*2);
dUl_dLM(xPosLinearIndLeftLeafI) = dUl_dLM(xPosLinearIndLeftLeafI) - ((leftLeafPosI-edges_l(xPosIndLeftLeafI)').^2)./((leftLeafPosM-leftLeafPosI).^2.*(widths(xPosIndLeftLeafI)').*2);
dUl_dLI(xPosLinearIndLeftLeafM) = dUl_dLI(xPosLinearIndLeftLeafM) - ((edges_r(xPosIndLeftLeafM)'-leftLeafPosM).^2)./((leftLeafPosM-leftLeafPosI).^2.*(widths(xPosIndLeftLeafM)').*2);
dUl_dLM(xPosLinearIndLeftLeafM) = dUl_dLM(xPosLinearIndLeftLeafM) + ((edges_r(xPosIndLeftLeafM)'-leftLeafPosM).*(leftLeafPosM+edges_r(xPosIndLeftLeafM)'-2*leftLeafPosI))./((leftLeafPosM-leftLeafPosI).^2.*(widths(xPosIndLeftLeafM)').*2);

dCr_dRI(xPosLinearIndRightLeafI) = dCr_dRI(xPosLinearIndRightLeafI) + ((rightLeafPosI-edges_l(xPosIndRightLeafI)').*(2*rightLeafPosM-rightLeafPosI-edges_l(xPosIndRightLeafI)'))./((rightLeafPosM-rightLeafPosI).^2.*(widths(xPosIndRightLeafI)').*2);
dCr_dRM(xPosLinearIndRightLeafI) = dCr_dRM(xPosLinearIndRightLeafI) - ((rightLeafPosI-edges_l(xPosIndRightLeafI)').^2)./((rightLeafPosM-rightLeafPosI).^2.*(widths(xPosIndRightLeafI)').*2);
dCr_dRI(xPosLinearIndRightLeafM) = dCr_dRI(xPosLinearIndRightLeafM) - ((edges_r(xPosIndRightLeafM)'-rightLeafPosM).^2)./((rightLeafPosM-rightLeafPosI).^2.*(widths(xPosIndRightLeafM)').*2);
dCr_dRM(xPosLinearIndRightLeafM) = dCr_dRM(xPosLinearIndRightLeafM) + ((edges_r(xPosIndRightLeafM)'-rightLeafPosM).*(rightLeafPosM+edges_r(xPosIndRightLeafM)'-2*rightLeafPosI))./((rightLeafPosM-rightLeafPosI).^2.*(widths(xPosIndRightLeafM)').*2);

for k = 1:n
    dUl_dLI(k,1:(xPosIndLeftLeafI(k)-1)) = 0;
    dUl_dLM(k,1:(xPosIndLeftLeafI(k)-1)) = 0;
    dUl_dLI(k,(xPosIndLeftLeafM(k)+1):numBix) = 0;
    dUl_dLM(k,(xPosIndLeftLeafM(k)+1):numBix) = 0;
    
    if xPosIndLeftLeafI(k) >= xPosIndLeftLeafM(k)
        % in discrete aperture, the xPosIndLeftLeafI is greater than
        % xPosIndLeftLeafM when leaf positions are at a bixel boundary
        
        %19 July 2017 in journal
        dUl_dLI(k,xPosIndLeftLeafI(k)) = -1/(2*widths(xPosIndLeftLeafI(k))');
        dUl_dLM(k,xPosIndLeftLeafM(k)) = -1/(2*widths(xPosIndLeftLeafM(k))');
        if leftLeafPosM(k)-leftLeafPosI(k) <= eps(max(lim_r))
            uncoveredByLeftLeaf(k,xPosIndLeftLeafI(k)) = (edges_r(xPosIndLeftLeafI(k))-leftLeafPosI(k))./widths(xPosIndLeftLeafI(k));
            uncoveredByLeftLeaf(k,xPosIndLeftLeafM(k)) = (edges_r(xPosIndLeftLeafM(k))-leftLeafPosM(k))./widths(xPosIndLeftLeafM(k));
        end
    end
    
    dCr_dRI(k,1:(xPosIndRightLeafI(k)-1)) = 0;
    dCr_dRM(k,1:(xPosIndRightLeafI(k)-1)) = 0;
    dCr_dRI(k,(xPosIndRightLeafM(k)+1):numBix) = 0;
    dCr_dRM(k,(xPosIndRightLeafM(k)+1):numBix) = 0;
    
    if xPosIndRightLeafI(k) >= xPosIndRightLeafM(k)
        dCr_dRI(k,xPosIndRightLeafI(k)) = -1/(2*widths(xPosIndRightLeafI(k))');
        dCr_dRM(k,xPosIndRightLeafM(k)) = -1/(2*widths(xPosIndRightLeafM(k))');
        if rightLeafPosM(k)-rightLeafPosI(k) <= eps(max(lim_r))
            coveredByRightLeaf(k,xPosIndRightLeafI(k)) = (edges_r(xPosIndRightLeafI(k))-rightLeafPosI(k))./widths(xPosIndRightLeafI(k));
            coveredByRightLeaf(k,xPosIndRightLeafM(k)) = (edges_r(xPosIndRightLeafM(k))-rightLeafPosM(k))./widths(xPosIndRightLeafM(k));
        end
    end
end

% store information for Jacobi preconditioning
sumGradSq = sumGradSq+(weightFactor_I).^2.*mean([sum((dUl_dLI).^2,2); sum((dUl_dLM.*weightFactor_F).^2,2); sum((dUl_dLM.*weightFactor_I).^2,2); sum((dCr_dRI).^2,2); sum((dCr_dRM.*weightFactor_F).^2,2); sum((dCr_dRM.*weightFactor_I).^2,2)]);

%% save the bixel weights
%fluence is equal to fluence not covered by left leaf minus
%fluence covered by left leaf
shapeMap = uncoveredByLeftLeaf-coveredByRightLeaf;
shapeMap = round2(shapeMap,15);
shapeMap(isnan(shapeMap)) = 0;

% find open bixels
%shapeMapIx = shapeMap > 0;
shapeMapIx = ~isnan(bixelIndMap);

currBixelIx = bixelIndMap(shapeMapIx);
w(currBixelIx) = w(currBixelIx) + shapeMap(shapeMapIx).*weight_I.*weightFactor_I;
shapeMapW = shapeMapW+shapeMap.*weight_I.*weightFactor_I;

%% save the gradients

if calcOptions.saveJacobian
    
    bixelJApVec_offset_temp = bixelJApVec_offset{phase_I};
    numSaveBixel = nnz(shapeMapIx);
    
    if calcOptions.DAOBeam
        % indices
        vectorIxMat_LI = repmat(vectorIx_LI',1,numBix);
        vectorIxMat_LF = repmat(vectorIx_LF',1,numBix);
        vectorIxMat_RI = repmat(vectorIx_RI',1,numBix);
        vectorIxMat_RF = repmat(vectorIx_RF',1,numBix);
        
        % wrt weight
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = shapeMap(shapeMapIx).*weightFactor_I./jacobiScale_I;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = DAOindex+(phase_I-1)*totalNumOfShapes;
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        
        % leftLeafPosM   = weightFactor_F.*leftLeafPosI+weightFactor_I.*leftLeafPosF;
        % rightLeafPosM  = weightFactor_F.*rightLeafPosI+weightFactor_I.*rightLeafPosF;
        
        % wrt initial left
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = dUl_dLI(shapeMapIx).*weight_I.*weightFactor_I;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LI(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        % extra weightFactor_F for I->M
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = dUl_dLM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_F;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LI(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        
        % wrt final left
        % extra weightFactor_I for F->M
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = dUl_dLM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_I;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LF(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        
        % wrt initial right
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = -dCr_dRI(shapeMapIx).*weight_I.*weightFactor_I;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RI(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        % extra weightFactor_F for I->M
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = -dCr_dRM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_F;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RI(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        
        % wrt final right
        % extra weightFactor_I for F->M
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = -dCr_dRM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_I;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RF(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        
    else
        % indices
        vectorIxMat_LF_last = repmat(vectorIx_LF_last',1,numBix);
        vectorIxMat_LI_next = repmat(vectorIx_LI_next',1,numBix);
        vectorIxMat_RF_last = repmat(vectorIx_RF_last',1,numBix);
        vectorIxMat_RI_next = repmat(vectorIx_RI_next',1,numBix);
        
        % leaf interpolation fractions/weights
        fracFromLastOptI = repmat(fracFromLastOptI,1,numBix);
        fracFromLastOptF = repmat(fracFromLastOptF,1,numBix);
        fracFromNextOptI = repmat(fracFromNextOptI,1,numBix);
        fracFromNextOptF = repmat(fracFromNextOptF,1,numBix);
        
        % wrt last weight
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromLastOpt*(time./time_last)*shapeMap(shapeMapIx).*weightFactor_I./jacobiScale_last_I;
        %bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.beam(i).doseAngleBordersDiff*fracFromLastOpt*updatedInfo.beam(apertureInfo.beam(i).lastOptIndex).gantryRot ...
        %/(updatedInfo.beam(apertureInfo.beam(i).lastOptIndex).doseAngleBordersDiff*updatedInfo.beam(i).gantryRot))*updatedInfo.beam(i).shape(j).shapeMap(shapeMapIx) ...
        %./ apertureInfo.beam(apertureInfo.beam(i).lastOptIndex).shape(1).jacobiScale;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = DAOindex_last+(phase_I-1)*totalNumOfShapes;
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        
        % wrt next weight
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = (1-fracFromLastOpt)*(time./time_next)*shapeMap(shapeMapIx).*weightFactor_I./jacobiScale_next_I;
        %bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.beam(i).doseAngleBordersDiff*(1-fracFromLastOpt)*updatedInfo.beam(apertureInfo.beam(i).nextOptIndex).gantryRot ...
        %/(updatedInfo.beam(apertureInfo.beam(i).nextOptIndex).doseAngleBordersDiff*updatedInfo.beam(i).gantryRot))*updatedInfo.beam(i).shape(j).shapeMap(shapeMapIx) ...
        %./ apertureInfo.beam(apertureInfo.beam(i).nextOptIndex).shape(1).jacobiScale;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = DAOindex_next+(phase_I-1)*totalNumOfShapes;
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        
        
        % updatedInfo.beam(i).shape(j).leftLeafPos_I = fracFromLastOptI.*leftLeafPos_F_last+fracFromNextOptI.*leftLeafPos_I_next;
        %updatedInfo.beam(i).shape(j).rightLeafPos_I = fracFromLastOptI.*rightLeafPos_F_last+fracFromNextOptI.*rightLeafPos_I_next;
        
        % updatedInfo.beam(i).shape(j).leftLeafPos_F = fracFromLastOptF.*leftLeafPos_F_last+fracFromNextOptF.*leftLeafPos_I_next;
        % updatedInfo.beam(i).shape(j).rightLeafPos_F = fracFromLastOptF.*rightLeafPos_F_last+fracFromNextOptF.*rightLeafPos_I_next;
        
        % leftLeafPosM   = weightFactor_F.*leftLeafPosI+weightFactor_I.*leftLeafPosF;
        % rightLeafPosM  = weightFactor_F.*rightLeafPosI+weightFactor_I.*rightLeafPosF;
        
        % wrt initial left (optimization vector)
        % initial (interpolated arc)
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromLastOptI(shapeMapIx).*dUl_dLI(shapeMapIx).*weight_I.*weightFactor_I;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LF_last(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        % extra weightFactor_F for I->M
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromLastOptI(shapeMapIx).*dUl_dLM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_F;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LF_last(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        % final (interpolated arc)
        % extra weightFactor_I for F->M
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromLastOptF(shapeMapIx).*dUl_dLM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_I;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LF_last(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        
        % wrt final left (optimization vector)
        % initial (interpolated arc)
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromNextOptI(shapeMapIx).*dUl_dLI(shapeMapIx).*weight_I.*weightFactor_I;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LI_next(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        % extra weightFactor_F for I->M
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromNextOptI(shapeMapIx).*dUl_dLM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_F;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LI_next(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        % final (interpolated arc)
        % extra weightFactor_I for F->M
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromNextOptF(shapeMapIx).*dUl_dLM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_I;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LI_next(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        
        % wrt initial right (optimization vector)
        % initial (interpolated arc)
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromLastOptI(shapeMapIx).*dCr_dRI(shapeMapIx).*weight_I.*weightFactor_I;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RF_last(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        % extra weightFactor_F for I->M
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromLastOptI(shapeMapIx).*dCr_dRM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_F;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RF_last(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        % final (interpolated arc)
        % extra weightFactor_I for F->M
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromLastOptF(shapeMapIx).*dCr_dRM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_I;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RF_last(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        
        % wrt final right (optimization vector)
        % initial (interpolated arc)
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromNextOptI(shapeMapIx).*dCr_dRI(shapeMapIx).*weight_I.*weightFactor_I;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RI_next(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        % extra weightFactor_F for I->M
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromNextOptI(shapeMapIx).*dCr_dRM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_F;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RI_next(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        % final (interpolated arc)
        % extra weightFactor_I for F->M
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromNextOptF(shapeMapIx).*dCr_dRM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_I;
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RI_next(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        
        % wrt last time
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = weightFactor_I.*(doseAngleBordersDiff.*timeFacCurr_last) ...
            .*(-fracFromLastDAO.*timeFracFromNextDAO.*(weight_last_I./doseAngleBordersDiff_next).*(time_next./time_last.^2) ...
            +(1-fracFromLastDAO).*timeFracFromLastDAO.*(weight_next_I./doseAngleBordersDiff_last).*(1./time_next)) ...
            * shapeMap(shapeMapIx);
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = tIx_last;
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
        
        % wrt next time
        bixelJApVec_vec(bixelJApVec_offset_temp+(1:numSaveBixel)) = weightFactor_I.*(doseAngleBordersDiff.*timeFacCurr_next) ...
            .*(fracFromLastDAO.*timeFracFromNextDAO.*(weight_last_I./doseAngleBordersDiff_next).*(1./time_last) ...
            -(1-fracFromLastDAO).*timeFracFromLastDAO.*(weight_next_I./doseAngleBordersDiff_last).*(time_last./time_next.^2)) ...
            * shapeMap(shapeMapIx);
        bixelJApVec_i(bixelJApVec_offset_temp+(1:numSaveBixel)) = tIx_next;
        bixelJApVec_j(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    end
    
    % update counters
    bixelJApVec_offset = bixelJApVec_offset_temp;
    counters.bixelJApVec_offset = bixelJApVec_offset;
    
end

end