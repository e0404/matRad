function jacobStruct = matRad_getJacobianStructure(optiProb,apertureInfoVec,dij,cst)
% matRad IPOPT callback: get jacobian structure for direct aperture optimization
% 
% call
%   jacobStruct = matRad_daoGetJacobStruct(optiProb,apertureInfoVec,dij,cst)
%
% input
%   optiProb:           option struct defining the type of optimization
%   apertureInfoVect:   aperture weights and shapes parameterized as vector
%   dij:                dose influence matrix
%   cst:                matRad cst struct
%
% output
%   jacobStruct: jacobian of constraint function
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

apertureInfo = optiProb.apertureInfo;

% jacobian structure of the dao constraints
% row indices
i = repmat(1:apertureInfo.totalNumOfLeafPairs,1,2);
% column indices
j = [apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs ...
     apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1:apertureInfo.totalNumOfShapes+2*apertureInfo.totalNumOfLeafPairs];

% -1 for left leaves, 1 for right leaves
s = ones(1,2*apertureInfo.totalNumOfLeafPairs);

jacobStruct_dao = sparse(i,j,s, ...
    apertureInfo.totalNumOfLeafPairs, ...
    apertureInfo.totalNumOfShapes+2*apertureInfo.totalNumOfLeafPairs, ...
    2*apertureInfo.totalNumOfLeafPairs);

jacobStruct_dos_bixel = matRad_getJacobianStructure@matRad_OptimizationProblem(optiProb,apertureInfo.bixelWeights,dij,cst);
% --> gives me a matrix with number of rows = num of constraints and tells
% me in th columns if a beamlet has an influence on this constraint

% for apertures I need to check if the very beam orientation of the aperture has a bixel
% that potentially influences the constraint

% for leaves I need to check if that particular leaf row has bixels that
% potentially influence the objective which works via apertureInfo.beam(i).bixelIndMap

% all stuff can be done per beam direction and then I use repmat to build
% up the big matrix

offset = 0;

if ~isempty(jacobStruct_dos_bixel)
    numOfConstraints = size(jacobStruct_dos_bixel,1);

    i_sparse = 1:numOfConstraints;
    i_sparse = kron(i_sparse,ones(1,numel(apertureInfo.apertureVector)));

    j_sparse = 1:numel(apertureInfo.apertureVector);
    j_sparse = repmat(j_sparse,1,numOfConstraints);

    jacobStructSparseVec = zeros(numOfConstraints*numel(apertureInfo.apertureVector),1);
    
    %counter = apertureInfo.totalNumOfShapes;
    for i = 1:numel(apertureInfo.beam)
    %for i = 1:size(apertureInfo.beam,2)
        
        % get used bixels in beam
        ixWeight = ~isnan(apertureInfo.beam(i).bixelIndMap);

        for j = 1:apertureInfo.beam(i).numOfShapes
            % first weight
            jacobStructSparseVec(offset+j == j_sparse) = jacobStructSparseVec(offset+j == j_sparse)+sum(jacobStruct_dos_bixel(:,apertureInfo.beam(i).bixelIndMap(ixWeight)),2);
            
            % now leaf positions
            for k = 1:apertureInfo.beam(i).numOfActiveLeafPairs
                
                ixLeaf  = ~isnan(apertureInfo.beam(i).bixelIndMap(k,:));
                indInBixVec = apertureInfo.beam(i).bixelIndMap(k,ixLeaf);
                
                indInOptVec = apertureInfo.beam(i).shape(1).vectorOffset+k-1+[0 apertureInfo.totalNumOfLeafPairs];
                indInSparseVec = repmat(indInOptVec,1,numOfConstraints)...
                    +kron((0:numOfConstraints-1)*numel(apertureInfo.apertureVector),ones(1,2));
                
                jacobStructSparseVec(indInSparseVec) = jacobStructSparseVec(indInSparseVec)+repelem(sum(jacobStruct_dos_bixel(:,indInBixVec),2),2,1);
            end
            
            offset = offset+1;
        end
    end
    
    jacobStructSparseVec(jacobStructSparseVec ~= 0) = 1;    
    jacobStruct_dos = sparse(i_sparse,j_sparse,jacobStructSparseVec,numOfConstraints,numel(apertureInfo.apertureVector));
    
else
    jacobStruct_dos = sparse(0,0);
end

% concatenate
jacobStruct = [jacobStruct_dos; jacobStruct_dao];
