function jacob = matRad_constraintJacobian(optiProb,apertureInfoVec,dij,cst)
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

apertureInfo = optiProb.apertureInfo;

% update apertureInfo if necessary
if ~isequal(apertureInfoVec,apertureInfo.apertureVector)
    apertureInfo = optiProb.matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVec);
end

% jacobian of the dao constraints

% row indices
i = repmat(1:apertureInfo.totalNumOfLeafPairs,1,2);
% column indices
j = [(apertureInfo.totalNumOfShapes+1):(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs) ...
     ((apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs)+1):(apertureInfo.totalNumOfShapes+2*apertureInfo.totalNumOfLeafPairs)];
    
% -1 for left leaves, 1 for right leaves
s = [-1*ones(1,apertureInfo.totalNumOfLeafPairs) ones(1,apertureInfo.totalNumOfLeafPairs)];

jacob_dao = sparse(i,j,s, ...
    apertureInfo.totalNumOfLeafPairs, ...
    numel(apertureInfoVec), ...
    2*apertureInfo.totalNumOfLeafPairs);

% compute jacobian of dosimetric constrainst

% dosimetric jacobian in bixel space
jacob_dos_bixel = matRad_constraintJacobian@matRad_OptimizationProblem(optiProb,apertureInfo.bixelWeights,dij,cst);

if ~isempty(jacob_dos_bixel)
    
    numOfConstraints = size(jacob_dos_bixel,1);
    
    i_sparse = 1:numOfConstraints;
    i_sparse = kron(i_sparse,ones(1,numel(apertureInfoVec)));
    
    j_sparse = 1:numel(apertureInfoVec);
    j_sparse = repmat(j_sparse,1,numOfConstraints);
    
    jacobSparseVec = zeros(numOfConstraints*size(apertureInfoVec,1),1);
    
    
    % 1. calculate jacobian for aperture weights
    % loop over all beams
    conOffset = 0;
    for i = 1:numel(apertureInfo.beam)
        
        % get used bixels in beam
        ix = ~isnan(apertureInfo.beam(i).bixelIndMap);
        
        % loop over all shapes and add up the gradients x openingFrac for this shape
        for j = 1:apertureInfo.beam(i).numOfShapes
            
            jacobSparseVec(conOffset+j == j_sparse) = jacob_dos_bixel(:,apertureInfo.beam(i).bixelIndMap(ix)) ...
                * apertureInfo.beam(i).shape(j).shapeMap(ix)./apertureInfo.beam(i).shape(j).jacobiScale;
        end
        
        % increment offset
        conOffset = conOffset + apertureInfo.beam(i).numOfShapes;
        
    end
    
    % 2. find corresponding bixel to the leaf Positions and aperture
    % weights to calculate the jacobian
    indInSparseVec = repmat(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2,1,numOfConstraints) ...
        +kron((0:numOfConstraints-1)*numel(apertureInfoVec),ones(1,apertureInfo.totalNumOfLeafPairs*2));
    
    jacobSparseVec(indInSparseVec) = ...
        reshape(transpose(( ones(numOfConstraints,1) * apertureInfoVec(apertureInfo.mappingMx(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2,2))' ) ...
        .* jacob_dos_bixel(:,apertureInfo.bixelIndices(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)) ./ ...
        (ones(numOfConstraints,1) * (apertureInfo.bixelWidth.*apertureInfo.jacobiScale(apertureInfo.mappingMx(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2,2)))')),[],1);
    
    
    % correct the sign for the left leaf positions
    indInSparseVec = repmat(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs,1,numOfConstraints) ...
        +kron((0:numOfConstraints-1)*numel(apertureInfoVec),ones(1,apertureInfo.totalNumOfLeafPairs));
    
    jacobSparseVec(indInSparseVec) = -jacobSparseVec(indInSparseVec);
    
    jacob_dos = sparse(i_sparse,j_sparse,jacobSparseVec,numOfConstraints,numel(apertureInfoVec));
    
else
    jacob_dos = sparse(0,0);
end

% concatenate
jacob = [jacob_dos; jacob_dao];
