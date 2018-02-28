function jacobStruct = matRad_daoGetJacobStruct(apertureInfo,dij,cst,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: get jacobian structure for direct aperture optimization
%
% call
%   jacobStruct = matRad_daoGetJacobStruct(apertureInfo,dij,cst)
%
% input
%   apertureInfo: aperture info struct
%   dij:          dose influence matrix
%   cst:          matRad cst struct
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
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    numel(apertureInfo.apertureVector), ...
    2*apertureInfo.totalNumOfLeafPairs);

jacobStruct_dos_bixel = matRad_getJacobStruct(dij,cst,options);
% --> gives me a matrix with number of rows = num of constraints and tells
% me in th columns if a beamlet has an influence on this constraint

% for apertures I need to check if the very beam orientation of the aperture has a bixel
% that potentially influences the constraint

% for leaves I need to check if that particular leaf row has bixels that
% potentially influence the objective which works via apertureInfo.beam(i).bixelIndMap

% all stuff can be done per beam direction and then I use repmat to build
% up the big matrix

if ~isempty(jacobStruct_dos_bixel)
    
    numOfConstraints = size(jacobStruct_dos_bixel,1);
    
    i_sparse = 1:numOfConstraints;
    i_sparse = kron(i_sparse,ones(1,numel(apertureInfo.apertureVector)));
    
    j_sparse = 1:numel(apertureInfo.apertureVector);
    j_sparse = repmat(j_sparse,1,numOfConstraints);
    
    jacobStructSparseVec = zeros(numOfConstraints*numel(apertureInfo.apertureVector),1);
    
    if apertureInfo.runVMAT
        
        offset = 1;
        DAOBeams = find([apertureInfo.propVMAT.beam.DAOBeam]);
        
        for i = 1:numel(apertureInfo.beam);
            
            % get used bixels in beam
            ixWeight = ~isnan(apertureInfo.beam(i).bixelIndMap);
            
            if apertureInfo.propVMAT.beam(i).DAOBeam
                % DAO beam, dont't worry about adding since this is just
                % struct, i.e. we are only interested if the element is
                % non-zero
                
                % first weight
                jacobStructSparseVec(offset == j_sparse) = sum(jacobStruct_dos_bixel(:,apertureInfo.beam(i).bixelIndMap(ixWeight)),2);
                
                % now leaf positions
                for k = 1:apertureInfo.beam(i).numOfActiveLeafPairs
                    
                    ixLeaf  = ~isnan(apertureInfo.beam(i).bixelIndMap(k,:));
                    indInBixVec = apertureInfo.beam(i).bixelIndMap(k,ixLeaf);
                    
                    indInOptVec = apertureInfo.beam(i).shape(1).vectorOffset+k-1+[0 apertureInfo.totalNumOfLeafPairs];
                    indInSparseVec = repmat(indInOptVec,1,numOfConstraints)...
                        +kron((0:numOfConstraints-1)*numel(apertureInfo.apertureVector),ones(1,2));
                    
                    jacobStructSparseVec(indInSparseVec) = repmat(sum(jacobStruct_dos_bixel(:,indInBixVec),2),2,1);
                end
                
                offset = offset+1;
            else
                % not DAO beam, these may contain bixels which affect the
                % constraints which are influenced by DAO leaf pairs that
                % do not affect the constraints (unlikely to happen, but it
                % might)
                
                %first weight
                
                %give fraction of gradient to previous optimized beam
                lastDAOInd = find(DAOBeams == apertureInfo.propVMAT.beam(i).lastDAOIndex,1);
                jacobStructSparseVec(lastDAOInd == j_sparse) = sum(jacobStruct_dos_bixel(:,apertureInfo.beam(i).bixelIndMap(ixWeight)),2);
                %give the other fraction to next optimized beam
                nextDAOInd = find(DAOBeams == apertureInfo.propVMAT.beam(i).nextDAOIndex,1);
                jacobStructSparseVec(nextDAOInd == j_sparse) = sum(jacobStruct_dos_bixel(:,apertureInfo.beam(i).bixelIndMap(ixWeight)),2);
                
                %now leaf pos
                
                for k = 1:apertureInfo.beam(i).numOfActiveLeafPairs
                    
                    ixLeaf  = ~isnan(apertureInfo.beam(i).bixelIndMap(k,:));
                    indInBixVec = apertureInfo.beam(i).bixelIndMap(k,ixLeaf);
                    
                    %give fraction of gradient to previous optimized beam
                    indInOptVec = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape(1).vectorOffset+k-1+[0 apertureInfo.totalNumOfLeafPairs];
                    indInSparseVec = repmat(indInOptVec,1,numOfConstraints)...
                        +kron((0:numOfConstraints-1)*numel(apertureInfo.apertureVector),ones(1,2));
                    
                    jacobStructSparseVec(indInSparseVec) = repmat(sum(jacobStruct_dos_bixel(:,indInBixVec),2),2,1);
                    
                    
                    %give the other fraction to next optimized beam
                    indInOptVec = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape(1).vectorOffset+k-1+[0 apertureInfo.totalNumOfLeafPairs];
                    indInSparseVec = repmat(indInOptVec,1,numOfConstraints)...
                        +kron((0:numOfConstraints-1)*numel(apertureInfo.apertureVector),ones(1,2));
                    
                    jacobStructSparseVec(indInSparseVec) = repmat(sum(jacobStruct_dos_bixel(:,indInBixVec),2),2,1);
                end
                
                %now time
                
                %give fraction of gradient to previous optimized beam
                lastDAOIndTime = lastDAOInd+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
                jacobStructSparseVec(lastDAOIndTime == j_sparse) = sum(jacobStruct_dos_bixel(:,apertureInfo.beam(i).bixelIndMap(ixWeight)),2);
                
                %give the other fraction to next optimized beam
                nextDAOIndTime = nextDAOInd+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
                jacobStructSparseVec(nextDAOIndTime == j_sparse) = sum(jacobStruct_dos_bixel(:,apertureInfo.beam(i).bixelIndMap(ixWeight)),2);
            end
        end
    else
        
        % new, correct
        offset = 0;
        for i = 1:numel(apertureInfo.beam);
            
            % get used bixels in beam
            ixWeight = ~isnan(apertureInfo.beam(i).bixelIndMap);
            
            for j = 1:apertureInfo.beam(i).numOfShapes
                % first weight
                jacobStructSparseVec(offset+j == j_sparse) = sum(jacobStruct_dos_bixel(:,apertureInfo.beam(i).bixelIndMap(ixWeight)),2);
                
                % now leaf positions
                for k = 1:apertureInfo.beam(i).numOfActiveLeafPairs
                    
                    ixLeaf  = ~isnan(apertureInfo.beam(i).bixelIndMap(k,:));
                    indInBixVec = apertureInfo.beam(i).bixelIndMap(k,ixLeaf);
                    
                    indInOptVec = apertureInfo.beam(i).shape(1).vectorOffset+k-1+[0 apertureInfo.totalNumOfLeafPairs];
                    indInSparseVec = repmat(indInOptVec,1,numOfConstraints)...
                        +kron((0:numOfConstraints-1)*numel(apertureInfo.apertureVector),ones(1,2));
                    
                    jacobStructSparseVec(indInSparseVec) = repmat(sum(jacobStruct_dos_bixel(:,indInBixVec),2),2,1);
                end
            end
            
            % increment offset
            offset = offset + apertureInfo.beam(i).numOfShapes;
        end
    end
    
    jacobStructSparseVec(jacobStructSparseVec ~= 0) = 1;
    
    jacobStruct_dos = sparse(i_sparse,j_sparse,jacobStructSparseVec,numOfConstraints,numel(apertureInfo.apertureVector));
else
    jacobStruct_dos = sparse(0,0);
end


if ~apertureInfo.runVMAT
    % concatenate
    jacobStruct = [jacobStruct_dao; jacobStruct_dos];
else
    
    i = sort(repmat(1:(apertureInfo.totalNumOfShapes-1),1,2));
    j = sort(repmat(1:apertureInfo.totalNumOfShapes,1,2));
    j(1) = [];
    j(end) = [];
    
    % get index values for the jacobian
    % variable index
    timeInd = (1+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2-1);
    currentLeftLeafInd = (apertureInfo.totalNumOfShapes+1):(apertureInfo.totalNumOfShapes+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd));
    currentRightLeafInd = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1):(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd));
    nextLeftLeafInd = (apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+1):(apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd));
    nextRightLeafInd = (apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1):(apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd));
    leftTimeInd = kron(j,ones(1,apertureInfo.beam(1).numOfActiveLeafPairs))+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
    rightTimeInd = kron(j,ones(1,apertureInfo.beam(1).numOfActiveLeafPairs))+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
    %leftTimeInd = repelem(j,apertureInfo.beam(1).numOfActiveLeafPairs)+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
    %rightTimeInd = repelem(j,apertureInfo.beam(1).numOfActiveLeafPairs)+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
    % constraint index
    constraintInd = 1:2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd);
    
    % jacobian of the leafspeed constraint
    i = repmat((i'-1)*apertureInfo.beam(1).numOfActiveLeafPairs,1,apertureInfo.beam(1).numOfActiveLeafPairs)+repmat(1:apertureInfo.beam(1).numOfActiveLeafPairs,2*numel(timeInd),1);
    i = reshape([i' i'+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd)],1,[]);
    
    i = [repmat(constraintInd,1,2) i];
    j = [currentLeftLeafInd currentRightLeafInd nextLeftLeafInd nextRightLeafInd leftTimeInd rightTimeInd];
    % first do jacob wrt current leaf position (left, right), then next leaf
    % position (left, right), then time (left, right)
    
    s = ones(1,numel(j));
    
    jacobStruct_lfspd = sparse(i,j,s,2*apertureInfo.beam(1).numOfActiveLeafPairs*(apertureInfo.totalNumOfShapes-1),numel(apertureInfo.apertureVector),numel(s));
    
    
    % jacobian of the doserate constraint
    i = repmat(1:apertureInfo.totalNumOfShapes,1,2);
    j = [1:apertureInfo.totalNumOfShapes (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+1):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2)];
    % first do jacob wrt weights, then wrt times
    
    s = ones(1,2*(apertureInfo.totalNumOfShapes));
    
    jacobStruct_dosrt = sparse(i,j,s,apertureInfo.totalNumOfShapes,numel(apertureInfo.apertureVector),2*apertureInfo.totalNumOfShapes);
    
    % concatenate
    jacobStruct = [jacobStruct_dao; jacobStruct_lfspd; jacobStruct_dosrt; jacobStruct_dos];
end


