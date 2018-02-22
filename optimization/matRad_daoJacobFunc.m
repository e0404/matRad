function jacob = matRad_daoJacobFunc(apertureInfoVec,dij,cst,options)
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
%   options:         option struct defining the type of optimization
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

% read in the global apertureInfo and apertureVector variables
global matRad_global_apertureInfo;
% update apertureInfo from the global variable
apertureInfo = matRad_global_apertureInfo;

% update apertureInfo, bixel weight vector an mapping of leafes to bixels
if ~isequal(apertureInfoVec,apertureInfo.apertureVector)
    apertureInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVec);
    matRad_global_apertureInfo = apertureInfo;
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
jacob_dos_bixel = matRad_jacobFuncWrapper(apertureInfo.bixelWeights,dij,cst,options);

if ~isempty(jacob_dos_bixel)
    
    numOfConstraints = size(jacob_dos_bixel,1);
    
    i_sparse = 1:numOfConstraints;
    i_sparse = kron(i_sparse,ones(1,numel(apertureInfoVec)));
    
    j_sparse = 1:numel(apertureInfoVec);
    j_sparse = repmat(j_sparse,1,numOfConstraints);
    
    jacobSparseVec = zeros(numOfConstraints*size(apertureInfoVec,1),1);
    
    if apertureInfo.runVMAT
        
        % 1. calculate aperatureGrad
        % 2. find corresponding bixel to the leaf Positions and aperture
        % weights to calculate the gradient
        
        offset = 1;
        DAOBeams = find([apertureInfo.propVMAT.beam.DAOBeam]);
        
        for i = 1:numel(apertureInfo.beam);
            
            % get used bixels in beam
            ix = ~isnan(apertureInfo.beam(i).bixelIndMap);
            
            if apertureInfo.propVMAT.beam(i).DAOBeam
                %optimized beam, do regular gradient
                %must always add to existing gradient, since gradient comes
                %from optimized and interpolated beams
                jacobSparseVec(offset == j_sparse) = jacobSparseVec(offset == j_sparse) + ...
                    jacob_dos_bixel(:,apertureInfo.beam(i).bixelIndMap(ix)) ...
                    * apertureInfo.beam(i).shape(1).shapeMap(ix)./apertureInfo.beam(i).shape(1).jacobiScale;
                
                %gradient wrt leaf positions
                indInOptVec = apertureInfo.beam(i).shape(1).vectorOffset-1+[(1:apertureInfo.beam(i).numOfActiveLeafPairs) apertureInfo.totalNumOfLeafPairs+(1:apertureInfo.beam(i).numOfActiveLeafPairs)];
                indInSparseVec = repmat(indInOptVec,1,numOfConstraints)...
                    +kron((0:numOfConstraints-1)*numel(apertureInfoVec),ones(1,apertureInfo.beam(i).numOfActiveLeafPairs*2));
                indInBixVec = apertureInfo.beam(i).bixOffset-1+[(1:apertureInfo.beam(i).numOfActiveLeafPairs) apertureInfo.doseTotalNumOfLeafPairs+(1:apertureInfo.beam(i).numOfActiveLeafPairs)];
                
                jacobSparseVec(indInSparseVec) = jacobSparseVec(indInSparseVec) + ...
                    reshape(transpose(apertureInfo.beam(i).shape(1).weight ...
                    .* jacob_dos_bixel(:,apertureInfo.bixelIndices(indInBixVec)) ./ apertureInfo.bixelWidth),[],1);
                
                offset = offset+1;
            else
                %not optimized beam, aperture weight is interpolated between
                %previous and next optimized weights
                
                %give fraction of gradient to previous optimized beam
                
                %first weight
                lastDAOInd = find(DAOBeams == apertureInfo.propVMAT.beam(i).lastDAOIndex,1);
                jacobSparseVec(lastDAOInd == j_sparse) = jacobSparseVec(lastDAOInd == j_sparse) + ...
                    apertureInfo.propVMAT.beam(i).fracFromLastDAO.*(apertureInfo.beam(i).time./apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).time) ...
                    * jacob_dos_bixel(:,apertureInfo.beam(i).bixelIndMap(ix)) * apertureInfo.beam(i).shape(1).shapeMap(ix)./apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape(1).jacobiScale;
                
                %now leaf pos
                indInOptVec = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape(1).vectorOffset-1+[(1:apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).numOfActiveLeafPairs) apertureInfo.totalNumOfLeafPairs+(1:apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).numOfActiveLeafPairs)];
                indInSparseVec = repmat(indInOptVec,1,numOfConstraints)...
                    +kron((0:numOfConstraints-1)*numel(apertureInfoVec),ones(1,apertureInfo.beam(i).numOfActiveLeafPairs*2));
                indInBixVec = apertureInfo.beam(i).bixOffset-1+[(1:apertureInfo.beam(i).numOfActiveLeafPairs) apertureInfo.doseTotalNumOfLeafPairs+(1:apertureInfo.beam(i).numOfActiveLeafPairs)];
                
                jacobSparseVec(indInSparseVec) = jacobSparseVec(indInSparseVec) + ...
                    apertureInfo.propVMAT.beam(i).fracFromLastDAO*reshape(transpose(apertureInfo.beam(i).shape(1).weight ...
                    .* jacob_dos_bixel(:,apertureInfo.bixelIndices(indInBixVec)) ./ apertureInfo.bixelWidth),[],1);
                
                %now time
                lastDAOIndTime = lastDAOInd+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
                jacobSparseVec(lastDAOIndTime == j_sparse) = jacobSparseVec(lastDAOIndTime == j_sparse) + ...
                    (apertureInfo.propVMAT.beam(i).doseAngleBordersDiff.*apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).timeFacCurr) ...
                    .*(-apertureInfo.propVMAT.beam(i).fracFromLastDAO.*apertureInfo.propVMAT.beam(i).timeFracFromNextDAO.*(apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape(1).weight./apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).doseAngleBordersDiff).*(apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).time./apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).time.^2) ...
                    +(1-apertureInfo.propVMAT.beam(i).fracFromLastDAO).*apertureInfo.propVMAT.beam(i).timeFracFromLastDAO.*(apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape(1).weight./apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).doseAngleBordersDiff).*(1./apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).time)) ...
                    * jacob_dos_bixel(:,apertureInfo.beam(i).bixelIndMap(ix)) * apertureInfo.beam(i).shape(1).shapeMap(ix);
                
                
                %give the other fraction to next optimized beam
                
                %first weight
                nextDAOInd = find(DAOBeams == apertureInfo.propVMAT.beam(i).nextDAOIndex,1);
                jacobSparseVec(nextDAOInd == j_sparse) = jacobSparseVec(nextDAOInd == j_sparse) + ...
                    (1-apertureInfo.propVMAT.beam(i).fracFromLastDAO).*(apertureInfo.beam(i).time./apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).time) ...
                    * jacob_dos_bixel(:,apertureInfo.beam(i).bixelIndMap(ix)) * apertureInfo.beam(i).shape(1).shapeMap(ix)./apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape(1).jacobiScale;
                
                %now leaf pos
                indInOptVec = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape(1).vectorOffset-1+[(1:apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).numOfActiveLeafPairs) apertureInfo.totalNumOfLeafPairs+(1:apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).numOfActiveLeafPairs)];
                indInSparseVec = repmat(indInOptVec,1,numOfConstraints)...
                    +kron((0:numOfConstraints-1)*numel(apertureInfoVec),ones(1,apertureInfo.beam(i).numOfActiveLeafPairs*2));
                indInBixVec = apertureInfo.beam(i).bixOffset-1+[(1:apertureInfo.beam(i).numOfActiveLeafPairs) apertureInfo.doseTotalNumOfLeafPairs+(1:apertureInfo.beam(i).numOfActiveLeafPairs)];
                
                jacobSparseVec(indInSparseVec) = jacobSparseVec(indInSparseVec) + ...
                    (1-apertureInfo.propVMAT.beam(i).fracFromLastDAO)*reshape(transpose(apertureInfo.beam(i).shape(1).weight ...
                    .* jacob_dos_bixel(:,apertureInfo.bixelIndices(indInBixVec)) ./ apertureInfo.bixelWidth),[],1);
                
                %now time
                nextDAOIndTime = nextDAOInd+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
                jacobSparseVec(nextDAOIndTime == j_sparse) = jacobSparseVec(nextDAOIndTime == j_sparse) + ...
                    (apertureInfo.propVMAT.beam(i).doseAngleBordersDiff.*apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).timeFacCurr) ...
                    .*(apertureInfo.propVMAT.beam(i).fracFromLastDAO.*apertureInfo.propVMAT.beam(i).timeFracFromNextDAO.*(apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape(1).weight./apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).doseAngleBordersDiff).*(1./apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).time) ...
                    -(1-apertureInfo.propVMAT.beam(i).fracFromLastDAO).*apertureInfo.propVMAT.beam(i).timeFracFromLastDAO.*(apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape(1).weight./apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).doseAngleBordersDiff).*(apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).time./apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).time.^2)) ...
                    * jacob_dos_bixel(:,apertureInfo.beam(i).bixelIndMap(ix)) * apertureInfo.beam(i).shape(1).shapeMap(ix);
            end
        end
    else
        
        % 1. calculate jacobian for aperture weights
        % loop over all beams
        offset = 0;
        for i = 1:numel(apertureInfo.beam);
            
            % get used bixels in beam
            ix = ~isnan(apertureInfo.beam(i).bixelIndMap);
            
            % loop over all shapes and add up the gradients x openingFrac for this shape
            for j = 1:apertureInfo.beam(i).numOfShapes
                
                jacobSparseVec(offset+j == j_sparse) = jacob_dos_bixel(:,apertureInfo.beam(i).bixelIndMap(ix)) ...
                    * apertureInfo.beam(i).shape(j).shapeMap(ix)./apertureInfo.beam(i).shape(j).jacobiScale;
            end
            
            % increment offset
            offset = offset + apertureInfo.beam(i).numOfShapes;
            
        end
        
        % 2. find corresponding bixel to the leaf Positions and aperture
        % weights to calculate the jacobian
        indInSparseVec = repmat(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2,1,numOfConstraints) ...
            +kron((0:numOfConstraints-1)*numel(apertureInfoVec),ones(1,apertureInfo.totalNumOfLeafPairs*2));
        
        jacobSparseVec(indInSparseVec) = ...
            reshape(transpose(( ones(numOfConstraints,1) * apertureInfoVec(apertureInfo.mappingMx(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2,2))' ) ...
            .* jacob_dos_bixel(:,apertureInfo.bixelIndices(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)) ./ ...
            (ones(numOfConstraints,1) * (apertureInfo.bixelWidth.*apertureInfo.jacobiScale(apertureInfo.mappingMx(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2,2)))')),[],1);
        
    end
    
    % correct the sign for the left leaf positions
    indInSparseVec = repmat(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs,1,numOfConstraints) ...
        +kron((0:numOfConstraints-1)*numel(apertureInfoVec),ones(1,apertureInfo.totalNumOfLeafPairs));
    
    jacobSparseVec(indInSparseVec) = -jacobSparseVec(indInSparseVec);
    
    jacob_dos = sparse(i_sparse,j_sparse,jacobSparseVec,numOfConstraints,numel(apertureInfoVec));
else
    jacob_dos = sparse(0,0);
end

if ~apertureInfo.runVMAT
    % concatenate
    jacob = [jacob_dao; jacob_dos];
else
    
    % get index values for the jacobian
    % variable index
    % value of constraints for leaves
    leftLeafPos  = apertureInfoVec((1:apertureInfo.totalNumOfLeafPairs)+apertureInfo.totalNumOfShapes);
    rightLeafPos = apertureInfoVec(1+apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2);
    
    % values of time differences of optimized gantry angles
    DAOInd = [apertureInfo.propVMAT.beam.DAOBeam];
    timeDAOBorderAngles = apertureInfoVec((1+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2):end);
    
    i = sort(repmat(1:(apertureInfo.totalNumOfShapes-1),1,2));
    j = sort(repmat(1:apertureInfo.totalNumOfShapes,1,2));
    j(1) = [];
    j(end) = [];
    
    timeFac = [apertureInfo.propVMAT.beam(DAOInd).timeFac]';
    timeFac(timeFac == 0) = [];
    
    timeFacMatrix = sparse(i,j,timeFac,(apertureInfo.totalNumOfShapes-1),apertureInfo.totalNumOfShapes);
    timeBNOptAngles = timeFacMatrix*timeDAOBorderAngles;
    
    currentLeftLeafInd = (apertureInfo.totalNumOfShapes+1):(apertureInfo.totalNumOfShapes+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles));
    currentRightLeafInd = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1):(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles));
    nextLeftLeafInd = (apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+1):(apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles));
    nextRightLeafInd = (apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1):(apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles));
    leftTimeInd = kron(j,ones(1,apertureInfo.beam(1).numOfActiveLeafPairs))+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
    rightTimeInd = kron(j,ones(1,apertureInfo.beam(1).numOfActiveLeafPairs))+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
    %leftTimeInd = repelem(j,apertureInfo.beam(1).numOfActiveLeafPairs)+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
    %rightTimeInd = repelem(j,apertureInfo.beam(1).numOfActiveLeafPairs)+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
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
    
    j_lfspd_t = -reshape([kron(abs(diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)),ones(1,2)).*repmat(timeFac',apertureInfo.beam(1).numOfActiveLeafPairs,1) ...
        kron(abs(diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)),ones(1,2)).*repmat(timeFac',apertureInfo.beam(1).numOfActiveLeafPairs,1)]./ ...
        repmat(kron((timeBNOptAngles.^2)',ones(1,2)),apertureInfo.beam(1).numOfActiveLeafPairs,2),[],1);
    
    s = [j_lfspd_cur; j_lfspd_nxt; j_lfspd_t];
    
    jacob_lfspd = sparse(i,j,s,2*apertureInfo.beam(1).numOfActiveLeafPairs*(apertureInfo.totalNumOfShapes-1),numel(apertureInfoVec),numel(s));
    
    
    % jacobian of the doserate constraint
    % values of doserate (MU/sec) between optimized gantry angles
    weights = apertureInfoVec(1:apertureInfo.totalNumOfShapes)./apertureInfo.jacobiScale;
    timeFacCurr = [apertureInfo.propVMAT.beam(DAOInd).timeFacCurr]';
    timeOptDoseBorderAngles = timeDAOBorderAngles.*timeFacCurr;
    
    i = repmat(1:apertureInfo.totalNumOfShapes,1,2);
    j = [1:apertureInfo.totalNumOfShapes (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+1):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2)];
    % first do jacob wrt weights, then wrt times
    
    s = [apertureInfo.weightToMU./(timeOptDoseBorderAngles.*apertureInfo.jacobiScale); -apertureInfo.weightToMU.*weights.*timeFacCurr./(timeOptDoseBorderAngles.^2)];
    
    jacob_dosrt = sparse(i,j,s,apertureInfo.totalNumOfShapes,numel(apertureInfoVec),2*(apertureInfo.totalNumOfShapes));
    
    % concatenate
    jacob = [jacob_dao; jacob_lfspd; jacob_dosrt; jacob_dos];
end


