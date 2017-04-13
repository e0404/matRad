function jacob = matRad_daoJacobFunc(apertureInfoVec,apertureInfo,dij,cst,type)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: jacobian function for direct aperture optimization
% 
% call
%   jacob = matRad_daoJacobFunc(apertureInfoVec,apertureInfo,dij,cst,type)
%
% input
%   apertureInfoVec: aperture info vector
%   apertureInfo:    aperture info struct
%   dij:             dose influence matrix
%   cst:             matRad cst struct
%   options: option struct defining the type of optimization
%
% output
%   jacob:           jacobian of constraint function
%
% References
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

% update apertureInfo if necessary
if ~isequal(apertureInfoVec,apertureInfo.apertureVector)
    apertureInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVec);
end

% jacobian of the dao constraints

% row indices
i = repmat(1:apertureInfo.totalNumOfLeafPairs,1,2);
% column indices
j = [apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs ...
     apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1:apertureInfo.totalNumOfShapes+2*apertureInfo.totalNumOfLeafPairs];

% -1 for left leaves, 1 for right leaves
s = [-1*ones(1,apertureInfo.totalNumOfLeafPairs) ones(1,apertureInfo.totalNumOfLeafPairs)];

jacob_dao = sparse(i,j,s, ...
    apertureInfo.totalNumOfLeafPairs, ...
    numel(apertureInfoVec), ...
    2*apertureInfo.totalNumOfLeafPairs);

% compute jacobian of dosimetric constrainst

% dosimetric jacobian in bixel space
jacob_dos_bixel = matRad_jacobFuncWrapper(apertureInfo.bixelWeights,dij,cst,type);

% allocate sparse matrix for dosimetric jacobian
jacob_dos = sparse(size(jacob_dos_bixel,1),numel(apertureInfoVec));

if ~isempty(jacob_dos)
    
    % 1. calculate jacobian for aperture weights
    % loop over all beams
    offset = 0;
    for i = 1:numel(apertureInfo.beam);

        % get used bixels in beam
        ix = ~isnan(apertureInfo.beam(i).bixelIndMap);

        % loop over all shapes and add up the gradients x openingFrac for this shape
        for j = 1:apertureInfo.beam(i).numOfShapes            
            jacob_dos(:,offset+j) = jacob_dos_bixel(:,apertureInfo.beam(i).bixelIndMap(ix)) ...
                                      * apertureInfo.beam(i).shape(j).shapeMap(ix);
        end

        % increment offset
        offset = offset + apertureInfo.beam(i).numOfShapes;

    end

    % 2. find corresponding bixel to the leaf Positions and aperture 
    % weights to calculate the jacobian
    jacob_dos(:,apertureInfo.totalNumOfShapes+1:end) = ...
            ( ones(size(jacob_dos,1),1) * apertureInfoVec(apertureInfo.mappingMx(apertureInfo.totalNumOfShapes+1:end,2))' ) ...
         .* jacob_dos_bixel(:,apertureInfo.bixelIndices(apertureInfo.totalNumOfShapes+1:end)) / apertureInfo.bixelWidth;

    % correct the sign for the left leaf positions
    jacob_dos(:,apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs) = ...
   -jacob_dos(:,apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs);

end

if ~isfield(apertureInfo.beam(1),'optimizeBeam')
    % concatenate
    jacob = [jacob_dao; jacob_dos];
else
    % get index values for the jacobian
    % variable index
    % value of constraints for leaves
    leftLeafPos  = apertureInfoVec([1:apertureInfo.totalNumOfLeafPairs]+apertureInfo.totalNumOfShapes);
    rightLeafPos = apertureInfoVec(1+apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2);
    
    % values of time differences of optimized gantry angles
    optInd = [apertureInfo.beam.optimizeBeam];
    timeOptBorderAngles = apertureInfoVec((1+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2):end);
    
    i = sort(repmat(1:(apertureInfo.totalNumOfShapes-1),1,2));
    j = sort(repmat(1:apertureInfo.totalNumOfShapes,1,2));
    j(1) = [];
    j(end) = [];
    
    timeFac = [apertureInfo.beam(optInd).timeFac]';
    timeFac(timeFac == 0) = [];
    
    timeFacMatrix = sparse(i,j,timeFac,(apertureInfo.totalNumOfShapes-1),apertureInfo.totalNumOfShapes);
    timeBNOptAngles = timeFacMatrix*timeOptBorderAngles;
    
    currentLeftLeafInd = (apertureInfo.totalNumOfShapes+1):(apertureInfo.totalNumOfShapes+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles));
    currentRightLeafInd = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1):(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles));
    nextLeftLeafInd = (apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+1):(apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles));
    nextRightLeafInd = (apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1):(apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles));
    leftTimeInd = repelem(j,apertureInfo.beam(1).numOfActiveLeafPairs)+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
    rightTimeInd = repelem(j,apertureInfo.beam(1).numOfActiveLeafPairs)+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
    % constraint index
    constraintInd = 1:2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles);
    
    
    % jacobian of the leafspeed constraint
    i = repmat((i'-1)*apertureInfo.beam(1).numOfActiveLeafPairs,1,apertureInfo.beam(1).numOfActiveLeafPairs)+repmat(1:apertureInfo.beam(1).numOfActiveLeafPairs,2*numel(timeBNOptAngles),1);
    i = reshape([i' i'+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles)],1,[]);
    
    i = [repmat(constraintInd,1,2) i];
    j = [currentLeftLeafInd currentRightLeafInd nextLeftLeafInd nextRightLeafInd leftTimeInd rightTimeInd];
    % first do jacob wrt current leaf position (left, right), then next leaf
    % position (left, right), then time (left, right)
    j_lfspd_cur = -reshape([sign(diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)) ...
        sign(diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2))]./ ...
        repmat(timeBNOptAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles),1);
    
    j_lfspd_nxt = reshape([sign(diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)) ...
        sign(diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2))]./ ...
        repmat(timeBNOptAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles),1);
    
    j_lfspd_t = -reshape([repelem(abs(diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)),1,2).*repmat(timeFac',apertureInfo.beam(1).numOfActiveLeafPairs,1) ...
        repelem(abs(diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)),1,2).*repmat(timeFac',apertureInfo.beam(1).numOfActiveLeafPairs,1)]./ ...
        repmat(repelem((timeBNOptAngles.^2)',2),apertureInfo.beam(1).numOfActiveLeafPairs,2),[],1);
        
    s = [j_lfspd_cur; j_lfspd_nxt; j_lfspd_t];
    
    jacob_lfspd = sparse(i,j,s,2*apertureInfo.beam(1).numOfActiveLeafPairs*(apertureInfo.totalNumOfShapes-1),numel(apertureInfoVec),numel(s));
    
    
    % jacobian of the doserate constraint
    % values of doserate (MU/sec) between optimized gantry angles
    weights = apertureInfoVec(1:apertureInfo.totalNumOfShapes);
    timeFacCurr = [apertureInfo.beam(optInd).timeFacCurr]';
    timeOptDoseBorderAngles = timeOptBorderAngles.*timeFacCurr;
    
    i = repmat(1:apertureInfo.totalNumOfShapes,1,2);
    j = [1:apertureInfo.totalNumOfShapes (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+1):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2)];
    % first do jacob wrt weights, then wrt times
    
    s = [apertureInfo.weightToMU.*timeOptDoseBorderAngles; -apertureInfo.weightToMU.*weights.*timeFacCurr./(timeOptDoseBorderAngles.^2)];
    
    jacob_dosrt = sparse(i,j,s,apertureInfo.totalNumOfShapes,numel(apertureInfoVec),2*(apertureInfo.totalNumOfShapes));
    
    % concatenate
    jacob = [jacob_dao; jacob_lfspd; jacob_dosrt; jacob_dos];
    
    %{
    IandFleftLeafPos  = apertureInfo.IandFapertureVector([1:apertureInfo.IandFtotalNumOfLeafPairs]+apertureInfo.totalNumOfShapes);
    IandFrightLeafPos = apertureInfo.IandFapertureVector(1+apertureInfo.IandFtotalNumOfLeafPairs+apertureInfo.totalNumOfShapes:apertureInfo.totalNumOfShapes+apertureInfo.IandFtotalNumOfLeafPairs*2);
    
    optInd = [apertureInfo.beam.optimizeBeam];
    timeOptBorderAngles = apertureInfoVec((1+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2):end);
    
    currIndGive = logical([apertureInfo.beam(optInd).timeFacCurr]);
    prevIndGive = logical([apertureInfo.beam(optInd).timeFacPrev]);
    nextIndGive = logical([apertureInfo.beam(optInd).timeFacNext]);
    
    timeFacCurr = [apertureInfo.beam(optInd).timeFacCurr]';
    timeFacCurr(~currIndGive) = [];
    timeFacPrev = [apertureInfo.beam(optInd).timeFacPrev]';
    timeFacPrev(~prevIndGive) = [];
    timeFacNext = [apertureInfo.beam(optInd).timeFacNext]';
    timeFacNext(~nextIndGive) = [];
    
    %IandFFacPrev = [apertureInfo.beam(optInd).IandFFacPrev]';
    %IandFFacPrev(~prevIndGive) = [];
    %IandFFacNext = [apertureInfo.beam(optInd).IandFFacNext]';
    %IandFFacNext(~nextIndGive) = [];
    
    timeInd = [apertureInfo.beam(optInd).IandFTimeInd];
    currInd = timeInd(2:3:end);
    prevInd = currInd(2:end)-1;
    prevInd(prevInd == currInd(1:(end-1))) = [];
    nextInd = currInd(1:(end-1))+1;
    nextInd(nextInd == currInd(2:end)) = [];
    
    timeIandFDoseBorderAngles = zeros(apertureInfo.numIandFBeam-1,1);
    timeIandFDoseBorderAngles(currInd) = timeIandFDoseBorderAngles(currInd) + timeOptBorderAngles(currIndGive).*timeFacCurr;
    timeIandFDoseBorderAngles(prevInd) = timeIandFDoseBorderAngles(prevInd) + timeOptBorderAngles(prevIndGive).*timeFacPrev;
    timeIandFDoseBorderAngles(nextInd) = timeIandFDoseBorderAngles(nextInd) + timeOptBorderAngles(nextIndGive).*timeFacNext;
    
    % values of average leaf speeds of optimized gantry angles
    c_lfspd = [abs(diff(reshape(IandFleftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2)) ...
        abs(diff(reshape(IandFrightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2))]./ ...
        repmat(timeIandFDoseBorderAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2);
    
    c_lfspd_sgn = [sign(diff(reshape(IandFleftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2)) ...
        sign(diff(reshape(IandFrightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2))]./ ...
        repmat(timeIandFDoseBorderAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2);
    
    jacob_lfspd_vec = zeros(6*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeIandFDoseBorderAngles),1);
    indInCon_vec = zeros(6*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeIandFDoseBorderAngles),1);
    indInVar_vec = zeros(6*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeIandFDoseBorderAngles),1);
    counter = 1;
    
    optCounter = 1;
    
    for i = 1:numel(apertureInfo.beam)
        if apertureInfo.beam(i).optimizeBeam
            
            indInIandF = apertureInfo.beam(i).IandFTimeInd;
            
            %leaf positions
            indInVar = reshape(repmat(apertureInfo.beam(i).shape(1).vectorOffset-1+[0 apertureInfo.totalNumOfLeafPairs],apertureInfo.beam(1).numOfActiveLeafPairs,1)+repmat((1:apertureInfo.beam(1).numOfActiveLeafPairs)',1,2),[],1);
            
            fac = apertureInfo.beam(i).IandFFac;
            
            if indInIandF(1) ~= 0
                indInC = [indInIandF(1) indInIandF(1)+numel(timeIandFDoseBorderAngles)];
                indInCon = reshape(repmat((indInC-1)*apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.beam(1).numOfActiveLeafPairs,1)+repmat((1:apertureInfo.beam(1).numOfActiveLeafPairs)',1,numel(indInC)),[],1);
                
                j = reshape(c_lfspd_sgn(:,indInC)*(fac(2)-fac(1)),[],1);
                
                indInVar_vec(counter:(counter+numel(j)-1)) = indInVar;
                indInCon_vec(counter:(counter+numel(j)-1)) = indInCon;
                jacob_lfspd_vec(counter:(counter+numel(j)-1)) = j+jacob_lfspd_vec(counter:(counter+numel(j)-1));
                counter = counter+numel(j);
            end
            
            indInC = [indInIandF(2) indInIandF(2)+numel(timeIandFDoseBorderAngles)];
            indInCon = reshape(repmat((indInC-1)*apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.beam(1).numOfActiveLeafPairs,1)+repmat((1:apertureInfo.beam(1).numOfActiveLeafPairs)',1,numel(indInC)),[],1);
            
            j = reshape(c_lfspd_sgn(:,indInC)*(fac(3)-fac(2)),[],1);
            
            indInVar_vec(counter:(counter+numel(j)-1)) = indInVar;
            indInCon_vec(counter:(counter+numel(j)-1)) = indInCon;
            jacob_lfspd_vec(counter:(counter+numel(j)-1)) = j+jacob_lfspd_vec(counter:(counter+numel(j)-1));
            counter = counter+numel(j);
            
            if indInIandF(3) ~= 0
                indInC = [indInIandF(3) indInIandF(3)+numel(timeIandFDoseBorderAngles)];
                indInCon = reshape(repmat((indInC-1)*apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.beam(1).numOfActiveLeafPairs,1)+repmat((1:apertureInfo.beam(1).numOfActiveLeafPairs)',1,numel(indInC)),[],1);
                
                j = reshape(c_lfspd_sgn(:,indInC)*(fac(4)-fac(3)),[],1);
                
                indInVar_vec(counter:(counter+numel(j)-1)) = indInVar;
                indInCon_vec(counter:(counter+numel(j)-1)) = indInCon;
                jacob_lfspd_vec(counter:(counter+numel(j)-1)) = j+jacob_lfspd_vec(counter:(counter+numel(j)-1));
                counter = counter+numel(j);
            end
            
            
            %times
            indInIandF(indInIandF == 0) = [];
            indInC = [indInIandF indInIandF+numel(timeIandFDoseBorderAngles)];
            
            indInVar = apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+optCounter;
            optCounter = optCounter+1;
            
            indInCon = reshape(repmat((indInC-1)*apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.beam(1).numOfActiveLeafPairs,1)+repmat((1:apertureInfo.beam(1).numOfActiveLeafPairs)',1,numel(indInC)),[],1);
            
            fac = [apertureInfo.beam(i).timeFacPrev apertureInfo.beam(i).timeFacCurr apertureInfo.beam(i).timeFacNext];
            fac(fac == 0) = [];
            
            j = -reshape(c_lfspd(:,indInC).*repmat(fac./timeIandFDoseBorderAngles(indInIandF)',apertureInfo.beam(1).numOfActiveLeafPairs,2),[],1);
            
            indInVar_vec(counter:(counter+numel(j)-1)) = indInVar;
            indInCon_vec(counter:(counter+numel(j)-1)) = indInCon;
            jacob_lfspd_vec(counter:(counter+numel(j)-1)) = j+jacob_lfspd_vec(counter:(counter+numel(j)-1));
            counter = counter+numel(j);
        end
    end
    
    jacob_lfspd_vec(indInVar_vec == 0) = [];
    indInVar_vec(indInVar_vec == 0) = [];
    indInCon_vec(indInCon_vec == 0) = [];
    
    jacob_lfspd = sparse(indInCon_vec,indInVar_vec,jacob_lfspd_vec,2*(apertureInfo.IandFtotalNumOfLeafPairs-apertureInfo.beam(1).numOfActiveLeafPairs),numel(apertureInfoVec),6*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeIandFDoseBorderAngles));
    %{
    currentLeftLeafInd = (apertureInfo.totalNumOfShapes+1):(apertureInfo.totalNumOfShapes+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeOptBorderAngles));
    currentRightLeafInd = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1):(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeOptBorderAngles));
    nextLeftLeafInd = (apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+1):(apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeOptBorderAngles));
    nextRightLeafInd = (apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1):(apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeOptBorderAngles));
    leftTimeInd = repelem((apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+1):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2),apertureInfo.beam(1).numOfActiveLeafPairs);
    rightTimeInd = repelem((apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+1):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2),apertureInfo.beam(1).numOfActiveLeafPairs);
    % constraint index
    constraintInd = 1:2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeIandFDoseBorderAngles);
    

    
    % jacobian of the leafspeed constraint
    i = repmat(constraintInd,1,3);
    j = [currentLeftLeafInd currentRightLeafInd nextLeftLeafInd nextRightLeafInd leftTimeInd rightTimeInd];
    % first do jacob wrt current leaf position (left, right), then next leaf
    % position (left, right), then time (left, right)
    
    %%%% first do jacob wrt prev leaf positions, then next
    
    j_lfspd_cur = -reshape([sign(diff(reshape(IandFleftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2)) ...
        sign(diff(reshape(IandFrightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2))]./ ...
        repmat(timeIandFDoseBorderAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeIandFDoseBorderAngles),1);
    
    j_lfspd_nxt = reshape([sign(diff(reshape(IandFleftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2)) ...
        sign(diff(reshape(IandFrightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2))]./ ...
        repmat(timeIandFDoseBorderAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeIandFDoseBorderAngles),1);
    
    j_lfspd_t = -reshape([abs(diff(reshape(IandFleftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2)) ...
        abs(diff(reshape(IandFrightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2))]./ ...
        repmat((timeIandFDoseBorderAngles.^2)',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeIandFDoseBorderAngles),1);
        
    s = [j_lfspd_cur; j_lfspd_nxt; j_lfspd_t];
    
    
    
    
    jacob_lfspd = sparse(i,j,s,2*(apertureInfo.IandFtotalNumOfLeafPairs-apertureInfo.beam(1).numOfActiveLeafPairs),numel(apertureInfoVec),6*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeIandFDoseBorderAngles));
    %}
    
    
    % jacobian of the doserate constraint
    % values of doserate (MU/sec) between optimized gantry angles
    weights = apertureInfoVec(1:apertureInfo.totalNumOfShapes);
    timeOptDoseBorderAngles = timeOptBorderAngles.*timeFacCurr;
    
    i = repmat(1:apertureInfo.totalNumOfShapes,1,2);
    j = [1:apertureInfo.totalNumOfShapes (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+1):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2)];
    % first do jacob wrt weights, then wrt times
    
    s = [apertureInfo.weightToMU.*timeOptDoseBorderAngles; -apertureInfo.weightToMU.*weights.*timeFacCurr./(timeOptDoseBorderAngles.^2)];
    
    jacob_dosrt = sparse(i,j,s,apertureInfo.totalNumOfShapes,numel(apertureInfoVec),2*(apertureInfo.totalNumOfShapes));
    
    % concatenate
    jacob = [jacob_dao; jacob_lfspd; jacob_dosrt; jacob_dos];
    %}
end


