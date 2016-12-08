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
%   type:   type of optimizaiton; either 'none','effect' or 'RBExD'
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
    timeInd = (1+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2-1);
    % value of constraints for leaves
    leftLeafPos  = apertureInfoVec([1:apertureInfo.totalNumOfLeafPairs]+apertureInfo.totalNumOfShapes);
    rightLeafPos = apertureInfoVec(1+apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2);
    % values of time differences of optimized gantry angles, used for
    % leafspeed and doserate constraint jacobians
    c_rottime = apertureInfoVec(timeInd);
    
    currentLeftLeafInd = (apertureInfo.totalNumOfShapes+1):(apertureInfo.totalNumOfShapes+apertureInfo.beam(1).numOfActiveLeafPairs*numel(c_rottime));
    currentRightLeafInd = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1):(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+apertureInfo.beam(1).numOfActiveLeafPairs*numel(c_rottime));
    nextLeftLeafInd = (apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+1):(apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.beam(1).numOfActiveLeafPairs*numel(c_rottime));
    nextRightLeafInd = (apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1):(apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+apertureInfo.beam(1).numOfActiveLeafPairs*numel(c_rottime));
    leftTimeInd = repelem((apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+1):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2-1),apertureInfo.beam(1).numOfActiveLeafPairs);
    rightTimeInd = repelem((apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+1):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2-1),apertureInfo.beam(1).numOfActiveLeafPairs);
    % constraint index
    constraintInd = 1:2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(c_rottime);
    

    
    % jacobian of the leafspeed constraint
    i = repmat(constraintInd,1,3);
    j = [currentLeftLeafInd currentRightLeafInd nextLeftLeafInd nextRightLeafInd leftTimeInd rightTimeInd];
    
    % first do jacob wrt current leaf position (left, right), then next leaf
    % position (left, right), then time (left, right)
    j_lfspd_cur = -reshape([sign(diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)) ...
        sign(diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2))]./ ...
        repmat(c_rottime',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(c_rottime),1);
    
    j_lfspd_nxt = reshape([sign(diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)) ...
        sign(diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2))]./ ...
        repmat(c_rottime',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(c_rottime),1);
    
    j_lfspd_t = -reshape([abs(diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)) ...
        abs(diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2))]./ ...
        repmat((c_rottime.^2)',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(c_rottime),1);
    
    %s = [repelem(-1./c_rottime,apertureInfo.beam(1).numOfActiveLeafPairs); repelem(-1./c_rottime,apertureInfo.beam(1).numOfActiveLeafPairs); ...
    %    repelem(1./c_rottime,apertureInfo.beam(1).numOfActiveLeafPairs); repelem(1./c_rottime,apertureInfo.beam(1).numOfActiveLeafPairs); j_lfspd_t];
    
    s = [j_lfspd_cur; j_lfspd_nxt; j_lfspd_t];
    
    jacob_lfspd = sparse(i,j,s,2*apertureInfo.beam(1).numOfActiveLeafPairs*(apertureInfo.totalNumOfShapes-1),numel(apertureInfoVec),6*apertureInfo.beam(1).numOfActiveLeafPairs*numel(c_rottime));
    
    
    % jacobian of the doserate constraint
    % values of doserate (MU/sec) between optimized gantry angles
    weights = apertureInfoVec(1:apertureInfo.totalNumOfShapes);
    weights(end) = [];
    optInd = find([apertureInfo.beam.optimizeBeam]);
    optInd(end) = [];
    nextOptAngleDiff = [apertureInfo.beam(optInd).nextOptAngleDiff]';
    nextAngleDiff = [apertureInfo.beam(optInd).nextAngleDiff]';
    
    i = repmat(1:(apertureInfo.totalNumOfShapes-1),1,2);
    j = [1:(apertureInfo.totalNumOfShapes-1) (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+1):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2-1)];
    % first do jacob wrt weights, then wrt times
    
    s = [apertureInfo.weightToMU.*(nextOptAngleDiff./c_rottime)./nextAngleDiff; -apertureInfo.weightToMU.*weights.*(nextOptAngleDiff./(c_rottime.^2))./nextAngleDiff];
    %%%IOAN%IARNAISDNUASIDUN LOOK AT THIS SQUARED TERM
    
    
    jacob_dosrt = sparse(i,j,s,apertureInfo.totalNumOfShapes-1,numel(apertureInfoVec),2*(apertureInfo.totalNumOfShapes-1));
    
    
    % concatenate
    jacob = [jacob_dao; jacob_lfspd; jacob_dosrt; jacob_dos];
end


