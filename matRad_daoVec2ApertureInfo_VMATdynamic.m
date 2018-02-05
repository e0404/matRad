function updatedInfo = matRad_daoVec2ApertureInfo_VMATdynamic(apertureInfo,apertureInfoVect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to translate the vector representation of the aperture
% shape and weight into an aperture info struct. At the same time, the
% updated bixel weight vector w is computed and a vector listing the
% correspondence between leaf tips and bixel indices for gradient
% calculation
%
% call
%   updatedInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVect)
%
% input
%   apertureInfo:     aperture shape info struct
%   apertureInfoVect: aperture weights and shapes parameterized as vector
%   touchingFlag:     if this is one, clean up instances of leaf touching,
%                     otherwise, do not
%
% output
%   updatedInfo: updated aperture shape info struct according to apertureInfoVect
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

% function to update the apertureInfo struct after the each iteraton of the
% optimization

w = zeros(apertureInfo.totalNumOfBixels,1);

% initializing variables
updatedInfo = apertureInfo;

updatedInfo.apertureVector = apertureInfoVect;


shapeInd = 1;
j = 1; %only 1 shape per beam for VMAT

% helper function to cope with numerical instabilities through rounding
round2 = @(a,b) round(a*10^b)/10^b;

% Jacobian matrix to be used in the DAO gradient function
% this tells us the gradient of a particular bixel with respect to an
% element in the apertureVector (aperture weight or leaf position)
% store as a vector for now, convert to sparse matrix later
bixelJApVec_sz = updatedInfo.totalNumOfOptBixels*5+(updatedInfo.totalNumOfBixels-updatedInfo.totalNumOfOptBixels)*10;
bixelJApVec_vec = zeros(1,bixelJApVec_sz);
% For optimized beams: 5 = (1 from weights) + (2 from left leaf positions (I and F)) + (2 from
% right leaf positions (I and F))
% For interpolated beams: multiply this number times 2 (influenced by the
% one before and the one after)

% vector indices
bixelJApVec_i = zeros(1,bixelJApVec_sz);
% bixel indices
bixelJApVec_j = zeros(1,bixelJApVec_sz);
% offset
bixelJApVec_offset = 0;



if ~all([updatedInfo.beam.optimizeBeam])
    for i = 1:numel(updatedInfo.beam)
        if updatedInfo.beam(i).optimizeBeam
            % update the shape weight
            updatedInfo.beam(i).shape(j).weight = apertureInfoVect(shapeInd)./updatedInfo.beam(i).shape(j).jacobiScale;
            
            updatedInfo.beam(i).MU = updatedInfo.beam(i).shape(j).weight*updatedInfo.weightToMU;
            updatedInfo.beam(i).time = apertureInfoVect(updatedInfo.totalNumOfShapes+updatedInfo.totalNumOfLeafPairs*2+shapeInd);
            updatedInfo.beam(i).gantryRot = updatedInfo.beam(i).optAngleBordersDiff/updatedInfo.beam(i).time;
            updatedInfo.beam(i).MURate = updatedInfo.beam(i).MU*updatedInfo.beam(i).gantryRot/updatedInfo.beam(i).doseAngleBordersDiff;
            
            shapeInd = shapeInd+1;
        end
    end
end

shapeInd = 1;
% loop over all beams
for i = 1:numel(updatedInfo.beam)
    
    %posOfRightCornerPixel = apertureInfo.beam(i).posOfCornerBixel(1) + (size(apertureInfo.beam(i).bixelIndMap,2)-1)*apertureInfo.bixelWidth;
    
    % pre compute left and right bixel edges
    edges_l = updatedInfo.beam(i).posOfCornerBixel(1)...
        + ([1:size(apertureInfo.beam(i).bixelIndMap,2)]-1-1/2)*updatedInfo.bixelWidth;
    edges_r = updatedInfo.beam(i).posOfCornerBixel(1)...
        + ([1:size(apertureInfo.beam(i).bixelIndMap,2)]-1+1/2)*updatedInfo.bixelWidth;
    
    centres = (edges_l+edges_r)/2;
    widths = edges_r-edges_l;
    numBix = size(apertureInfo.beam(i).bixelIndMap,2);
    
    % get dimensions of 2d matrices that store shape/bixel information
    n = apertureInfo.beam(i).numOfActiveLeafPairs;
    
    % extract leaf position and aperture weight information from vector
    % for optimized beams, this is straightforward
    % for non-optimized beams, interpolate MURate and leaf speed
    if updatedInfo.beam(i).optimizeBeam
        % update the shape weight
        updatedInfo.beam(i).shape(j).weight = apertureInfoVect(shapeInd)./updatedInfo.beam(i).shape(j).jacobiScale;
        
        updatedInfo.beam(i).MU = updatedInfo.beam(i).shape(j).weight*updatedInfo.weightToMU;
        updatedInfo.beam(i).time = apertureInfoVect(updatedInfo.totalNumOfShapes+updatedInfo.totalNumOfLeafPairs*2+shapeInd);
        updatedInfo.beam(i).gantryRot = updatedInfo.beam(i).optAngleBordersDiff/updatedInfo.beam(i).time;
        updatedInfo.beam(i).MURate = updatedInfo.beam(i).MU*updatedInfo.beam(i).gantryRot/updatedInfo.beam(i).doseAngleBordersDiff;
        
        
        vectorIx_LI = updatedInfo.beam(i).shape(j).vectorOffset(1) + ([1:n]-1);
        vectorIx_RI = vectorIx_LI+apertureInfo.totalNumOfLeafPairs;
        leftLeafPos_I = apertureInfoVect(vectorIx_LI);
        rightLeafPos_I = apertureInfoVect(vectorIx_RI);
        
        vectorIx_LF = updatedInfo.beam(i).shape(j).vectorOffset(2) + ([1:n]-1);
        vectorIx_RF = vectorIx_LF+apertureInfo.totalNumOfLeafPairs;
        leftLeafPos_F = apertureInfoVect(vectorIx_LF);
        rightLeafPos_F = apertureInfoVect(vectorIx_RF);
        
        % update information in shape structure
        updatedInfo.beam(i).shape(j).leftLeafPos_I  = leftLeafPos_I;
        updatedInfo.beam(i).shape(j).rightLeafPos_I = rightLeafPos_I;
        
        updatedInfo.beam(i).shape(j).leftLeafPos_F  = leftLeafPos_F;
        updatedInfo.beam(i).shape(j).rightLeafPos_F = rightLeafPos_F;
    else
        %MURate is interpolated between MURates of optimized apertures
        updatedInfo.beam(i).MURate = updatedInfo.beam(i).fracFromLastOpt*updatedInfo.beam(updatedInfo.beam(i).lastOptIndex).MURate+(1-updatedInfo.beam(i).fracFromLastOpt)*updatedInfo.beam(updatedInfo.beam(i).nextOptIndex).MURate;
        updatedInfo.beam(i).gantryRot = 1./(updatedInfo.beam(i).timeFracFromLastOpt./updatedInfo.beam(updatedInfo.beam(i).lastOptIndex).gantryRot+updatedInfo.beam(i).timeFracFromNextOpt./updatedInfo.beam(updatedInfo.beam(i).nextOptIndex).gantryRot);
        
        updatedInfo.beam(i).MU = updatedInfo.beam(i).MURate*updatedInfo.beam(i).doseAngleBordersDiff/updatedInfo.beam(i).gantryRot;
        updatedInfo.beam(i).shape(1).weight = updatedInfo.beam(i).MU./updatedInfo.weightToMU;
        
        vectorIx_LI_last = updatedInfo.beam(updatedInfo.beam(i).lastOptIndex).shape(j).vectorOffset(1) + ([1:n]-1);
        vectorIx_RI_last = vectorIx_LI_last+apertureInfo.totalNumOfLeafPairs;
        leftLeafPos_I_last = apertureInfoVect(vectorIx_LI_last);
        rightLeafPos_I_last = apertureInfoVect(vectorIx_RI_last);
        
        vectorIx_LI_next = updatedInfo.beam(updatedInfo.beam(i).nextOptIndex).shape(j).vectorOffset(1) + ([1:n]-1);
        vectorIx_RI_next = vectorIx_LI_next+apertureInfo.totalNumOfLeafPairs;
        leftLeafPos_I_next = apertureInfoVect(vectorIx_LI_next);
        rightLeafPos_I_next = apertureInfoVect(vectorIx_RI_next);
        
        vectorIx_LF_last = updatedInfo.beam(updatedInfo.beam(i).lastOptIndex).shape(j).vectorOffset(2) + ([1:n]-1);
        vectorIx_RF_last = vectorIx_LF_last+apertureInfo.totalNumOfLeafPairs;
        leftLeafPos_F_last = apertureInfoVect(vectorIx_LF_last);
        rightLeafPos_F_last = apertureInfoVect(vectorIx_RF_last);
        
        vectorIx_LF_next = updatedInfo.beam(updatedInfo.beam(i).nextOptIndex).shape(j).vectorOffset(2) + ([1:n]-1);
        vectorIx_RF_next = vectorIx_LF_next+apertureInfo.totalNumOfLeafPairs;
        leftLeafPos_F_next = apertureInfoVect(vectorIx_LF_next);
        rightLeafPos_F_next = apertureInfoVect(vectorIx_RF_next);
        
        
        updatedInfo.beam(i).shape(j).leftLeafPos_I = updatedInfo.beam(i).fracFromLastOptI*leftLeafPos_I_last+(1-updatedInfo.beam(i).fracFromLastOptI)*leftLeafPos_I_next;
        updatedInfo.beam(i).shape(j).rightLeafPos_I = updatedInfo.beam(i).fracFromLastOptI*rightLeafPos_I_last+(1-updatedInfo.beam(i).fracFromLastOptI)*rightLeafPos_I_next;
        
        updatedInfo.beam(i).shape(j).leftLeafPos_F = updatedInfo.beam(i).fracFromLastOptF*leftLeafPos_F_last+(1-updatedInfo.beam(i).fracFromLastOptF)*leftLeafPos_F_next;
        updatedInfo.beam(i).shape(j).rightLeafPos_F = updatedInfo.beam(i).fracFromLastOptF*rightLeafPos_F_last+(1-updatedInfo.beam(i).fracFromLastOptF)*rightLeafPos_F_next;
    end
    

    
    % set the initial leaf positions to the minimum leaf positions
    % always, instead of the leaf positions at the actual beginning
    % of the arc
    % this simplifies the calculation
    % remember which one is actually I and F in leftMinInd
    [leftLeafPosI,leftMinInd] = min([updatedInfo.beam(i).shape(j).leftLeafPos_I,updatedInfo.beam(i).shape(j).leftLeafPos_F],[],2);
    leftLeafPosF = max([updatedInfo.beam(i).shape(j).leftLeafPos_I,updatedInfo.beam(i).shape(j).leftLeafPos_F],[],2);
    [rightLeafPosI,rightMinInd] = min([updatedInfo.beam(i).shape(j).rightLeafPos_I,updatedInfo.beam(i).shape(j).rightLeafPos_F],[],2);
    rightLeafPosF = max([updatedInfo.beam(i).shape(j).rightLeafPos_I,updatedInfo.beam(i).shape(j).rightLeafPos_F],[],2);
    
    if updatedInfo.beam(i).optimizeBeam
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
        tempL = vectorIx_LI_last;
        tempR = vectorIx_RI_last;
        vectorIx_LI_last(leftMinInd == 2) = vectorIx_LF_last(leftMinInd == 2);
        vectorIx_LF_last(leftMinInd == 2) = tempL(leftMinInd == 2);
        vectorIx_RI_last(rightMinInd == 2) = vectorIx_RF_last(rightMinInd == 2);
        vectorIx_RF_last(rightMinInd == 2) = tempR(rightMinInd == 2);
        
        tempL = vectorIx_LI_next;
        tempR = vectorIx_RI_next;
        vectorIx_LI_next(leftMinInd == 2) = vectorIx_LF_next(leftMinInd == 2);
        vectorIx_LF_next(leftMinInd == 2) = tempL(leftMinInd == 2);
        vectorIx_RI_next(rightMinInd == 2) = vectorIx_RF_next(rightMinInd == 2);
        vectorIx_RF_next(rightMinInd == 2) = tempR(rightMinInd == 2);
    end
    
    %{
            leftLeafPosI = round2(leftLeafPosI,10);
            leftLeafPosF = round2(leftLeafPosF,10);
            rightLeafPosI = round2(rightLeafPosI,10);
            rightLeafPosF = round2(rightLeafPosF,10);
    %}
    leftLeafPosI(leftLeafPosI <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(leftLeafPosI <= apertureInfo.beam(i).lim_l);
    leftLeafPosF(leftLeafPosF <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(leftLeafPosF <= apertureInfo.beam(i).lim_l);
    rightLeafPosI(rightLeafPosI <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(rightLeafPosI <= apertureInfo.beam(i).lim_l);
    rightLeafPosF(rightLeafPosF <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(rightLeafPosF <= apertureInfo.beam(i).lim_l);
    leftLeafPosI(leftLeafPosI >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(leftLeafPosI >= apertureInfo.beam(i).lim_r);
    leftLeafPosF(leftLeafPosF >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(leftLeafPosF >= apertureInfo.beam(i).lim_r);
    rightLeafPosI(rightLeafPosI >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(rightLeafPosI >= apertureInfo.beam(i).lim_r);
    rightLeafPosF(rightLeafPosF >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(rightLeafPosF >= apertureInfo.beam(i).lim_r);
    
    % find bixel indices where leaves are located
    xPosIndLeftLeafI = min(floor((leftLeafPosI-edges_l(1))./apertureInfo.bixelWidth)+1,numBix);
    xPosIndLeftLeafF = max(ceil((leftLeafPosF-edges_r(1))./apertureInfo.bixelWidth)+1,1);
    xPosIndRightLeafI = min(floor((rightLeafPosI-edges_l(1))./apertureInfo.bixelWidth)+1,numBix);
    xPosIndRightLeafF = max(ceil((rightLeafPosF-edges_r(1))./apertureInfo.bixelWidth)+1,1);
    %
    xPosLinearIndLeftLeafI = sub2ind([n numBix],(1:n)',xPosIndLeftLeafI);
    xPosLinearIndLeftLeafF = sub2ind([n numBix],(1:n)',xPosIndLeftLeafF);
    xPosLinearIndRightLeafI = sub2ind([n numBix],(1:n)',xPosIndRightLeafI);
    xPosLinearIndRightLeafF = sub2ind([n numBix],(1:n)',xPosIndRightLeafF);
    
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
    
    
    %% gradients
    %21 June 2017
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
        if leftLeafPosI(k) == leftLeafPosF(k)
            %19 July 2017 in journal
            dUl_dLI(k,xPosIndLeftLeafI(k)) = -1/(2*widths(xPosIndLeftLeafI(k))');
            dUl_dLF(k,xPosIndLeftLeafF(k)) = -1/(2*widths(xPosIndLeftLeafF(k))');
            uncoveredByLeftLeaf(k,xPosIndLeftLeafI(k)) = (edges_r(xPosIndLeftLeafI(k))-leftLeafPosI(k))./widths(xPosIndLeftLeafI(k));
            uncoveredByLeftLeaf(k,xPosIndLeftLeafF(k)) = (edges_r(xPosIndLeftLeafF(k))-leftLeafPosF(k))./widths(xPosIndLeftLeafF(k));
        end
        
        dCr_dRI(k,1:(xPosIndRightLeafI(k)-1)) = 0;
        dCr_dRF(k,1:(xPosIndRightLeafI(k)-1)) = 0;
        dCr_dRI(k,(xPosIndRightLeafF(k)+1):numBix) = 0;
        dCr_dRF(k,(xPosIndRightLeafF(k)+1):numBix) = 0;
        if rightLeafPosI(k) == rightLeafPosF(k)
            dCr_dRI(k,xPosIndRightLeafI(k)) = -1/(2*widths(xPosIndRightLeafI(k))');
            dCr_dRF(k,xPosIndRightLeafF(k)) = -1/(2*widths(xPosIndRightLeafF(k))');
            coveredByRightLeaf(k,xPosIndRightLeafI(k)) = (edges_r(xPosIndRightLeafI(k))-rightLeafPosI(k))./widths(xPosIndRightLeafI(k));
            coveredByRightLeaf(k,xPosIndRightLeafF(k)) = (edges_r(xPosIndRightLeafF(k))-rightLeafPosF(k))./widths(xPosIndRightLeafF(k));
        end
    end
    
    %% save the bixel weights
    %fluence is equal to fluence not covered by left leaf minus
    %fluence covered by left leaf
    tempMap = uncoveredByLeftLeaf-coveredByRightLeaf;
    tempMap = round2(tempMap,15);
    tempMap(isnan(tempMap)) = 0;
    
    % find open bixels
    tempMapIx = tempMap > 0;
    
    currBixelIx = apertureInfo.beam(i).bixelIndMap(tempMapIx);
    w(currBixelIx) = w(currBixelIx) + tempMap(tempMapIx)*updatedInfo.beam(i).shape(j).weight;
    
    % save the tempMap (we need to apply a positivity operator !)
    updatedInfo.beam(i).shape(j).shapeMap = (tempMap  + abs(tempMap))  / 2;
    
    %% save the gradients
    
    if updatedInfo.beam(i).optimizeBeam
        % indices
        saveBixelIx = ~isnan(apertureInfo.beam(i).bixelIndMap);
        numSaveBixel = nnz(saveBixelIx);
        vectorIx_LI = repmat(vectorIx_LI',1,numBix);
        vectorIx_LF = repmat(vectorIx_LF',1,numBix);
        vectorIx_RI = repmat(vectorIx_RI',1,numBix);
        vectorIx_RF = repmat(vectorIx_RF',1,numBix);
        
        % wrt weight
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = updatedInfo.beam(i).shape(j).shapeMap(saveBixelIx)./apertureInfo.beam(i).shape(1).jacobiScale;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = shapeInd;
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt initial left
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = dUl_dLI(saveBixelIx)*updatedInfo.beam(i).shape(j).weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LI(saveBixelIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt final left
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = dUl_dLF(saveBixelIx)*updatedInfo.beam(i).shape(j).weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LF(saveBixelIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt initial right
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = -dCr_dRI(saveBixelIx)*updatedInfo.beam(i).shape(j).weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RI(saveBixelIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt final right
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = -dCr_dRF(saveBixelIx)*updatedInfo.beam(i).shape(j).weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RF(saveBixelIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        
        if apertureInfo.updateJacobi
            dijScaleFactor = mean(apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes)./apertureInfo.jacobiScale)/(2*apertureInfo.bixelWidth);
            
            if apertureInfo.jacobi
                % "incorrect"
                %updatedInfo.beam(i).shape(j).jacobiScale = sqrt(sum(updatedInfo.beam(i).shape(j).shapeMap(:)));
                % "correct"
                updatedInfo.beam(i).shape(j).jacobiScale = (dijScaleFactor./apertureInfo.beam(i).shape(j).weight).*sqrt(sum(updatedInfo.beam(i).shape(j).shapeMap(:).^2))./sqrt(mean([sum(dUl_dLI.^2,2); sum(dUl_dLF.^2,2); sum(dCr_dRI.^2,2); sum(dCr_dRF.^2,2)]));
            end
            % rescale the vector from the weight using the current
            % iteration scaling factor
            apertureInfoVect(shapeInd) = updatedInfo.beam(i).shape(j).jacobiScale*updatedInfo.beam(i).shape(j).weight;
            
            updatedInfo.jacobiScale(shapeInd) = updatedInfo.beam(i).shape(j).jacobiScale;
        end
        
        % increment shape index
        shapeInd = shapeInd +1;
    else
        % indices
        saveBixelIx = ~isnan(apertureInfo.beam(i).bixelIndMap);
        numSaveBixel = nnz(saveBixelIx);
        
        % first do last optimized beam
        vectorIx_LI_last = repmat(vectorIx_LI_last',1,numBix);
        vectorIx_LF_last = repmat(vectorIx_LF_last',1,numBix);
        vectorIx_RI_last = repmat(vectorIx_RI_last',1,numBix);
        vectorIx_RF_last = repmat(vectorIx_RF_last',1,numBix);
        
        % wrt weight
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.beam(i).doseAngleBordersDiff*updatedInfo.beam(i).fracFromLastOpt/updatedInfo.beam(apertureInfo.beam(i).lastOptIndex).doseAngleBordersDiff)*updatedInfo.beam(i).shape(j).shapeMap(saveBixelIx);
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = nnz([updatedInfo.beam(1:updatedInfo.beam(i).lastOptIndex).optimizeBeam]);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt initial left
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = updatedInfo.beam(i).fracFromLastOptI*dUl_dLI(saveBixelIx)*updatedInfo.beam(i).shape(j).weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LI_last(saveBixelIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt final left
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = updatedInfo.beam(i).fracFromLastOptF*dUl_dLF(saveBixelIx)*updatedInfo.beam(i).shape(j).weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LF_last(saveBixelIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt initial right
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = -updatedInfo.beam(i).fracFromLastOptI*dCr_dRI(saveBixelIx)*updatedInfo.beam(i).shape(j).weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RI_last(saveBixelIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt final right
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = -updatedInfo.beam(i).fracFromLastOptF*dCr_dRF(saveBixelIx)*updatedInfo.beam(i).shape(j).weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RF_last(saveBixelIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        
        % now do next optimized beam
        vectorIx_LI_next = repmat(vectorIx_LI_next',1,numBix);
        vectorIx_LF_next = repmat(vectorIx_LF_next',1,numBix);
        vectorIx_RI_next = repmat(vectorIx_RI_next',1,numBix);
        vectorIx_RF_next = repmat(vectorIx_RF_next',1,numBix);
        
        % wrt weight
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.beam(i).doseAngleBordersDiff*(1-updatedInfo.beam(i).fracFromLastOpt)/updatedInfo.beam(apertureInfo.beam(i).nextOptIndex).doseAngleBordersDiff)*updatedInfo.beam(i).shape(j).shapeMap(saveBixelIx);
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = nnz([updatedInfo.beam(1:updatedInfo.beam(i).nextOptIndex).optimizeBeam]);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt initial left
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (1-updatedInfo.beam(i).fracFromLastOptI)*dUl_dLI(saveBixelIx)*updatedInfo.beam(i).shape(j).weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LI_next(saveBixelIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt final left
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (1-updatedInfo.beam(i).fracFromLastOptF)*dUl_dLF(saveBixelIx)*updatedInfo.beam(i).shape(j).weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LF_next(saveBixelIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt initial right
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = -(1-updatedInfo.beam(i).fracFromLastOptI)*dCr_dRI(saveBixelIx)*updatedInfo.beam(i).shape(j).weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RI_next(saveBixelIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt final right
        bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = -(1-updatedInfo.beam(i).fracFromLastOptF)*dCr_dRF(saveBixelIx)*updatedInfo.beam(i).shape(j).weight;
        bixelJApVec_i(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RF_next(saveBixelIx);
        bixelJApVec_j(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
    end
    
end

updatedInfo.bixelWeights = w;
updatedInfo.apertureVector = apertureInfoVect;
updatedInfo.bixelJApVec = sparse(bixelJApVec_i,bixelJApVec_j,bixelJApVec_vec,numel(apertureInfoVect),updatedInfo.totalNumOfBixels,bixelJApVec_sz);

end