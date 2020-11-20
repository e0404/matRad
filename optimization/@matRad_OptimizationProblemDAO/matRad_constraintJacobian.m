function jacob = matRad_constraintJacobian(optiProb,apertureInfoVec,dij,cst)
% matRad IPOPT callback: jacobian function for direct aperture optimization
% 
% call
%   jacob = matRad_daoJacobFunc(optiProb,apertureInfoVec,dij,cst)
%
% input
%   optiProb:        option struct defining the type of optimization
%   apertureInfoVec: aperture info vector
%   dij:             dose influence matrix
%   cst:             matRad cst struct
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
j = [apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs ...
     apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1:apertureInfo.totalNumOfShapes+2*apertureInfo.totalNumOfLeafPairs];

% -1 for left leaves, 1 for right leaves
s = [-1*ones(1,apertureInfo.totalNumOfLeafPairs) ones(1,apertureInfo.totalNumOfLeafPairs)];

jacob_dao = sparse(i,j,s, ...
    apertureInfo.totalNumOfLeafPairs, ...
    apertureInfo.totalNumOfShapes+2*apertureInfo.totalNumOfLeafPairs, ...
    2*apertureInfo.totalNumOfLeafPairs);

% compute jacobian of dosimetric constrainst

% dosimetric jacobian in bixel space
jacob_dos_bixel = matRad_constraintJacobian@matRad_OptimizationProblem(optiProb,apertureInfo.bixelWeights,dij,cst);

% allocate sparse matrix for dosimetric jacobian
jacob_dos = sparse(size(jacob_dos_bixel,1),numel(apertureInfoVec));

if ~isempty(jacob_dos)
    
    % 1. calculate jacobian for aperture weights
    % loop over all beams
    offset = 0;
    for i = 1:numel(apertureInfo.beam)

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

% concatenate
jacob = [jacob_dao;jacob_dos];
