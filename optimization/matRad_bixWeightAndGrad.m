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

weight        = variables.weight;
leftLeafPos_I   = variables.leftLeafPos_I;
leftLeafPos_F   = variables.leftLeafPos_F;
rightLeafPos_I  = variables.rightLeafPos_I;
rightLeafPos_F  = variables.rightLeafPos_F;

bixelJApVec_offset = counters.bixelJApVec_offset;


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
    
    if calcOptions.DAOBeam
        jacobiScale = variables.jacobiScale;
        
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
        
        weight_last = variables.weight_last;
        weight_next = variables.weight_next;
        jacobiScale_last = variables.jacobiScale_last;
        jacobiScale_next = variables.jacobiScale_next;
        
        time_last = variables.time_last;
        time_next = variables.time_next;
        time = variables.time;
        
        fracFromLastOptI = variables.fracFromLastOptI;
        fracFromLastOptF = variables.fracFromLastOptF;
        fracFromNextOptI = variables.fracFromNextOptI;
        fracFromNextOptF = variables.fracFromNextOptF;
        fracFromLastOpt = variables.fracFromLastOpt;
        
        % replicate
        fracFromLastOptI = repmat(fracFromLastOptI,1,numBix);
        fracFromLastOptF = repmat(fracFromLastOptF,1,numBix);
        fracFromNextOptI = repmat(fracFromNextOptI,1,numBix);
        fracFromNextOptF = repmat(fracFromNextOptF,1,numBix);
        
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

% find bixel indices where leaves are located
xPosIndLeftLeafI = min(floor((leftLeafPosI-edges_l(1))./bixelWidth)+1,numBix);
xPosIndLeftLeafF = min(floor((leftLeafPosF-edges_l(1))./bixelWidth)+1,numBix);
xPosIndRightLeafI = min(floor((rightLeafPosI-edges_l(1))./bixelWidth)+1,numBix);
xPosIndRightLeafF = min(floor((rightLeafPosF-edges_l(1))./bixelWidth)+1,numBix);
%
xPosLinearIndLeftLeafI = sub2ind([n numBix],(1:n)',xPosIndLeftLeafI);
xPosLinearIndLeftLeafF = sub2ind([n numBix],(1:n)',xPosIndLeftLeafF);
xPosLinearIndRightLeafI = sub2ind([n numBix],(1:n)',xPosIndRightLeafI);
xPosLinearIndRightLeafF = sub2ind([n numBix],(1:n)',xPosIndRightLeafF);


%
% leaves sweep from _I to _F, with weight
%


%% bixel weight calculation

%calculate fraction of fluence uncovered by left leaf
%initial computation
uncoveredByLeftLeaf = bsxfun(@minus,centres,leftLeafPosI)./repmat(leftLeafPosF-leftLeafPosI,1,numBix);
%correct for overshoot in initial and final leaf positions
uncoveredByLeftLeaf(xPosLinearIndLeftLeafI) = uncoveredByLeftLeaf(xPosLinearIndLeftLeafI) + (leftLeafPosI-edges_l(xPosIndLeftLeafI)').^2./((leftLeafPosF-leftLeafPosI).*(widths(xPosIndLeftLeafI)').*2);
uncoveredByLeftLeaf(xPosLinearIndLeftLeafF) = uncoveredByLeftLeaf(xPosLinearIndLeftLeafF) - (edges_r(xPosIndLeftLeafF)'-leftLeafPosF).^2./((leftLeafPosF-leftLeafPosI).*(widths(xPosIndLeftLeafF)').*2);
%round <0 to 0, >1 to 1
uncoveredByLeftLeaf(uncoveredByLeftLeaf < 0) = 0;
uncoveredByLeftLeaf(uncoveredByLeftLeaf > 1) = 1;

%calculate fraction of fluence covered by right leaf
%initial computation
coveredByRightLeaf = bsxfun(@minus,centres,rightLeafPosI)./repmat(rightLeafPosF-rightLeafPosI,1,numBix);
%correct for overshoot in initial and final leaf positions
coveredByRightLeaf(xPosLinearIndRightLeafI) = coveredByRightLeaf(xPosLinearIndRightLeafI) + (rightLeafPosI-edges_l(xPosIndRightLeafI)').^2./((rightLeafPosF-rightLeafPosI).*(widths(xPosIndRightLeafI)').*2);
coveredByRightLeaf(xPosLinearIndRightLeafF) = coveredByRightLeaf(xPosLinearIndRightLeafF) - (edges_r(xPosIndRightLeafF)'-rightLeafPosF).^2./((rightLeafPosF-rightLeafPosI).*(widths(xPosIndRightLeafF)').*2);
%round <0 to 0, >1 to 1
coveredByRightLeaf(coveredByRightLeaf < 0) = 0;
coveredByRightLeaf(coveredByRightLeaf > 1) = 1;

%% gradient calculation

dUl_dLI = bsxfun(@minus,centres,leftLeafPosF)./(repmat(leftLeafPosF-leftLeafPosI,1,numBix)).^2;
dUl_dLF = bsxfun(@minus,leftLeafPosI,centres)./(repmat(leftLeafPosF-leftLeafPosI,1,numBix)).^2;

dCr_dRI = bsxfun(@minus,centres,rightLeafPosF)./(repmat(rightLeafPosF-rightLeafPosI,1,numBix)).^2;
dCr_dRF = bsxfun(@minus,rightLeafPosI,centres)./(repmat(rightLeafPosF-rightLeafPosI,1,numBix)).^2;

dUl_dLI(xPosLinearIndLeftLeafI) = dUl_dLI(xPosLinearIndLeftLeafI) + ((leftLeafPosI-edges_l(xPosIndLeftLeafI)').*(2*leftLeafPosF-leftLeafPosI-edges_l(xPosIndLeftLeafI)'))./((leftLeafPosF-leftLeafPosI).^2.*(widths(xPosIndLeftLeafI)').*2);
dUl_dLF(xPosLinearIndLeftLeafI) = dUl_dLF(xPosLinearIndLeftLeafI) - ((leftLeafPosI-edges_l(xPosIndLeftLeafI)').^2)./((leftLeafPosF-leftLeafPosI).^2.*(widths(xPosIndLeftLeafI)').*2);
dUl_dLI(xPosLinearIndLeftLeafF) = dUl_dLI(xPosLinearIndLeftLeafF) - ((edges_r(xPosIndLeftLeafF)'-leftLeafPosF).^2)./((leftLeafPosF-leftLeafPosI).^2.*(widths(xPosIndLeftLeafF)').*2);
dUl_dLF(xPosLinearIndLeftLeafF) = dUl_dLF(xPosLinearIndLeftLeafF) + ((edges_r(xPosIndLeftLeafF)'-leftLeafPosF).*(leftLeafPosF+edges_r(xPosIndLeftLeafF)'-2*leftLeafPosI))./((leftLeafPosF-leftLeafPosI).^2.*(widths(xPosIndLeftLeafF)').*2);

dCr_dRI(xPosLinearIndRightLeafI) = dCr_dRI(xPosLinearIndRightLeafI) + ((rightLeafPosI-edges_l(xPosIndRightLeafI)').*(2*rightLeafPosF-rightLeafPosI-edges_l(xPosIndRightLeafI)'))./((rightLeafPosF-rightLeafPosI).^2.*(widths(xPosIndRightLeafI)').*2);
dCr_dRF(xPosLinearIndRightLeafI) = dCr_dRF(xPosLinearIndRightLeafI) - ((rightLeafPosI-edges_l(xPosIndRightLeafI)').^2)./((rightLeafPosF-rightLeafPosI).^2.*(widths(xPosIndRightLeafI)').*2);
dCr_dRI(xPosLinearIndRightLeafF) = dCr_dRI(xPosLinearIndRightLeafF) - ((edges_r(xPosIndRightLeafF)'-rightLeafPosF).^2)./((rightLeafPosF-rightLeafPosI).^2.*(widths(xPosIndRightLeafF)').*2);
dCr_dRF(xPosLinearIndRightLeafF) = dCr_dRF(xPosLinearIndRightLeafF) + ((edges_r(xPosIndRightLeafF)'-rightLeafPosF).*(rightLeafPosF+edges_r(xPosIndRightLeafF)'-2*rightLeafPosI))./((rightLeafPosF-rightLeafPosI).^2.*(widths(xPosIndRightLeafF)').*2);

for k = 1:n
    dUl_dLI(k,1:(xPosIndLeftLeafI(k)-1)) = 0;
    dUl_dLF(k,1:(xPosIndLeftLeafI(k)-1)) = 0;
    dUl_dLI(k,(xPosIndLeftLeafF(k)+1):numBix) = 0;
    dUl_dLF(k,(xPosIndLeftLeafF(k)+1):numBix) = 0;
    
    if xPosIndLeftLeafI(k) >= xPosIndLeftLeafF(k)
        % in discrete aperture, the xPosIndLeftLeafI is greater than
        % xPosIndLeftLeafM when leaf positions are at a bixel boundary
        
        %19 July 2017 in journal
        dUl_dLI(k,xPosIndLeftLeafI(k)) = -1/(2*widths(xPosIndLeftLeafI(k))');
        dUl_dLF(k,xPosIndLeftLeafF(k)) = -1/(2*widths(xPosIndLeftLeafF(k))');
        if leftLeafPosF(k)-leftLeafPosI(k) <= eps(max(lim_r))
            uncoveredByLeftLeaf(k,xPosIndLeftLeafI(k)) = (edges_r(xPosIndLeftLeafI(k))-leftLeafPosI(k))./widths(xPosIndLeftLeafI(k));
            uncoveredByLeftLeaf(k,xPosIndLeftLeafF(k)) = (edges_r(xPosIndLeftLeafF(k))-leftLeafPosF(k))./widths(xPosIndLeftLeafF(k));
        end
    end
    
    dCr_dRI(k,1:(xPosIndRightLeafI(k)-1)) = 0;
    dCr_dRF(k,1:(xPosIndRightLeafI(k)-1)) = 0;
    dCr_dRI(k,(xPosIndRightLeafF(k)+1):numBix) = 0;
    dCr_dRF(k,(xPosIndRightLeafF(k)+1):numBix) = 0;
    
    if xPosIndRightLeafI(k) >= xPosIndRightLeafF(k)
        dCr_dRI(k,xPosIndRightLeafI(k)) = -1/(2*widths(xPosIndRightLeafI(k))');
        dCr_dRF(k,xPosIndRightLeafF(k)) = -1/(2*widths(xPosIndRightLeafF(k))');
        if rightLeafPosF(k)-rightLeafPosI(k) <= eps(max(lim_r))
            coveredByRightLeaf(k,xPosIndRightLeafI(k)) = (edges_r(xPosIndRightLeafI(k))-rightLeafPosI(k))./widths(xPosIndRightLeafI(k));
            coveredByRightLeaf(k,xPosIndRightLeafF(k)) = (edges_r(xPosIndRightLeafF(k))-rightLeafPosF(k))./widths(xPosIndRightLeafF(k));
        end
    end
end

% store information for Jacobi preconditioning
sumGradSq = sumGradSq+mean([sum((dUl_dLI).^2,2); sum((dUl_dLF).^2,2); sum((dUl_dLF).^2,2); sum((dCr_dRI).^2,2); sum((dCr_dRF).^2,2); sum((dCr_dRF).^2,2)]);

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
w(currBixelIx) = w(currBixelIx) + shapeMap(shapeMapIx).*weight;
shapeMapW = shapeMapW+shapeMap.*weight;

%% save the gradients

if calcOptions.saveJacobian
    
    numSaveBixel = nnz(shapeMapIx);
    
    if calcOptions.DAOBeam
        % indices
        vectorIxMat_LI = repmat(vectorIx_LI',1,numBix);
        vectorIxMat_LF = repmat(vectorIx_LF',1,numBix);
        vectorIxMat_RI = repmat(vectorIx_RI',1,numBix);
        vectorIxMat_RF = repmat(vectorIx_RF',1,numBix);
        
        % wrt weight
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = shapeMap(shapeMapIx)./jacobiScale;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = DAOindex;
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt initial left
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = dUl_dLI(shapeMapIx).*weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIxMat_LI(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt final left
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = dUl_dLF(shapeMapIx).*weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIxMat_LF(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt initial right
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = -dCr_dRI(shapeMapIx).*weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIxMat_RI(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt final right
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = -dCr_dRF(shapeMapIx).*weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIxMat_RF(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
    else
        % indices
        vectorIxMat_LF_last = repmat(vectorIx_LF_last',1,numBix);
        vectorIxMat_LI_next = repmat(vectorIx_LI_next',1,numBix);
        vectorIxMat_RF_last = repmat(vectorIx_RF_last',1,numBix);
        vectorIxMat_RI_next = repmat(vectorIx_RI_next',1,numBix);
        
        % wrt last weight
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = fracFromLastOpt*(time./time_last)*shapeMap(shapeMapIx)./jacobiScale_last;
        %bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.beam(i).doseAngleBordersDiff*fracFromLastOpt*updatedInfo.beam(apertureInfo.beam(i).lastOptIndex).gantryRot ...
        %/(updatedInfo.beam(apertureInfo.beam(i).lastOptIndex).doseAngleBordersDiff*updatedInfo.beam(i).gantryRot))*updatedInfo.beam(i).shape(j).shapeMap(shapeMapIx) ...
        %./ apertureInfo.beam(apertureInfo.beam(i).lastOptIndex).shape(1).jacobiScale;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = DAOindex_last;
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt next weight
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (1-fracFromLastOpt)*(time./time_next)*shapeMap(shapeMapIx)./jacobiScale_next;
        %bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.beam(i).doseAngleBordersDiff*(1-fracFromLastOpt)*updatedInfo.beam(apertureInfo.beam(i).nextOptIndex).gantryRot ...
        %/(updatedInfo.beam(apertureInfo.beam(i).nextOptIndex).doseAngleBordersDiff*updatedInfo.beam(i).gantryRot))*updatedInfo.beam(i).shape(j).shapeMap(shapeMapIx) ...
        %./ apertureInfo.beam(apertureInfo.beam(i).nextOptIndex).shape(1).jacobiScale;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = DAOindex_next;
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        
        % updatedInfo.beam(i).shape(j).leftLeafPos_I = fracFromLastOptI.*leftLeafPos_F_last+fracFromNextOptI.*leftLeafPos_I_next;
        %updatedInfo.beam(i).shape(j).rightLeafPos_I = fracFromLastOptI.*rightLeafPos_F_last+fracFromNextOptI.*rightLeafPos_I_next;
        
        % updatedInfo.beam(i).shape(j).leftLeafPos_F = fracFromLastOptF.*leftLeafPos_F_last+fracFromNextOptF.*leftLeafPos_I_next;
        % updatedInfo.beam(i).shape(j).rightLeafPos_F = fracFromLastOptF.*rightLeafPos_F_last+fracFromNextOptF.*rightLeafPos_I_next;
        
        % wrt initial left (optimization vector)
        % initial (interpolated arc)
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = fracFromLastOptI(shapeMapIx).*dUl_dLI(shapeMapIx).*weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIxMat_LF_last(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        % final (interpolated arc)
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = fracFromLastOptF(shapeMapIx).*dUl_dLF(shapeMapIx).*weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIxMat_LF_last(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt final left (optimization vector)
        % initial (interpolated arc)
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = fracFromNextOptI(shapeMapIx).*dUl_dLI(shapeMapIx).*weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIxMat_LI_next(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        % final (interpolated arc)
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = fracFromNextOptF(shapeMapIx).*dUl_dLF(shapeMapIx).*weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIxMat_LI_next(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt initial right (optimization vector)
        % initial (interpolated arc)
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = -fracFromLastOptI(shapeMapIx).*dCr_dRI(shapeMapIx).*weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIxMat_RF_last(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        % final (interpolated arc)
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = -fracFromLastOptF(shapeMapIx).*dCr_dRF(shapeMapIx).*weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIxMat_RF_last(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt final right (optimization vector)
        % initial (interpolated arc)
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = -fracFromNextOptI(shapeMapIx).*dCr_dRI(shapeMapIx).*weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIxMat_RI_next(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        % final (interpolated arc)
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = -fracFromNextOptF(shapeMapIx).*dCr_dRF(shapeMapIx).*weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIxMat_RI_next(shapeMapIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt last time
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (doseAngleBordersDiff.*timeFacCurr_last) ...
            .*(-fracFromLastDAO.*timeFracFromNextDAO.*(weight_last./doseAngleBordersDiff_next).*(time_next./time_last.^2) ...
            +(1-fracFromLastDAO).*timeFracFromLastDAO.*(weight_next./doseAngleBordersDiff_last).*(1./time_next)) ...
            * shapeMap(shapeMapIx);
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = tIx_last;
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt next time
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (doseAngleBordersDiff.*timeFacCurr_next) ...
            .*(fracFromLastDAO.*timeFracFromNextDAO.*(weight_last./doseAngleBordersDiff_next).*(1./time_last) ...
            -(1-fracFromLastDAO).*timeFracFromLastDAO.*(weight_next./doseAngleBordersDiff_last).*(time_last./time_next.^2)) ...
            * shapeMap(shapeMapIx);
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = tIx_next;
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
    end
    
end

% update counters
counters.bixelJApVec_offset = bixelJApVec_offset;

end