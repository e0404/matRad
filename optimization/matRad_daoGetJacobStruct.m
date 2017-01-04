function jacobStruct = matRad_daoGetJacobStruct(apertureInfo,dij,cst)
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
    apertureInfo.totalNumOfShapes+2*apertureInfo.totalNumOfLeafPairs, ...
    2*apertureInfo.totalNumOfLeafPairs);

jacobStruct_dos_bixel = matRad_getJacobStruct(dij,cst);
% --> gives me a matrix with number of rows = num of constraints and tells
% me in th columns if a beamlet has an influence on this constraint

% for apertures I need to check if the very beam orientation of the aperture has a bixel
% that potentially influences the constraint

% for leaves I need to check if that particular leaf row has bixels that
% potentially influence the objective which works via apertureInfo.beam(i).bixelIndMap

% all stuff can be done per beam direction and then I use repmat to build
% up the big matrix

% allocate
jacobStruct_dos = sparse(size(jacobStruct_dos_bixel,1),size(jacobStruct_dao,2));

if ~isempty(jacobStruct_dos)
    
    % first aperture weights
    for i = 1:apertureInfo.totalNumOfShapes
        currBeam             = apertureInfo.mappingMx(i,1);
        currBixelIxInBeam    = dij.beamNum == currBeam;
        jacobStruct_dos(:,i) = spones(sum(jacobStruct_dos_bixel(:,currBixelIxInBeam),2));
    end

    % second leaves
    counter = 0;
    for i = 1:size(apertureInfo.beam,2)
        for j = 1:apertureInfo.beam(i).numOfShapes
            for k = 1:apertureInfo.beam(i).numOfActiveLeafPairs
                counter = counter + 1;
                bixelIxInCurrRow = ~isnan(apertureInfo.beam(i).bixelIndMap(k,:));
                jacobStruct_dos(:,counter+[0 apertureInfo.totalNumOfLeafPairs]) = ...     
                    repmat(spones(sum(jacobStruct_dos_bixel(:,bixelIxInCurrRow),2)),1,2);
            end
        end    
    end

end



if ~isfield(apertureInfo.beam(1),'optimizeBeam')
    % concatenate
    jacobStruct = [jacobStruct_dao; jacobStruct_dos];
else
    % get index values for the jacobian
    % variable index
    timeInd = (1+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2-1);
    currentLeftLeafInd = (apertureInfo.totalNumOfShapes+1):(apertureInfo.totalNumOfShapes+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd));
    currentRightLeafInd = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1):(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd));
    nextLeftLeafInd = (apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+1):(apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd));
    nextRightLeafInd = (apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1):(apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd));
    leftTimeInd = repelem((apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+1):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2-1),apertureInfo.beam(1).numOfActiveLeafPairs);
    rightTimeInd = repelem((apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+1):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2-1),apertureInfo.beam(1).numOfActiveLeafPairs);
    % constraint index
    constraintInd = 1:2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd);
    
    % jacobian of the leafspeed constraint
    i = repmat(constraintInd,1,3);
    j = [currentLeftLeafInd currentRightLeafInd nextLeftLeafInd nextRightLeafInd leftTimeInd rightTimeInd];
    % first do jacob wrt current leaf position (left, right), then next leaf
    % position (left, right), then time (left, right)
    
    s = ones(1,6*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd));
    
    jacobStruct_lfspd = sparse(i,j,s,2*apertureInfo.beam(1).numOfActiveLeafPairs*(apertureInfo.totalNumOfShapes-1),numel(apertureInfo.apertureVector),6*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd));
    
    
    % jacobian of the doserate constraint
    i = repmat(1:(apertureInfo.totalNumOfShapes-1),1,2);
    j = [1:(apertureInfo.totalNumOfShapes-1) (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+1):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2-1)];
    % first do jacob wrt weights, then wrt times
    
    s = ones(1,2*(apertureInfo.totalNumOfShapes-1));
    
    jacobStruct_dosrt = sparse(i,j,s,apertureInfo.totalNumOfShapes-1,numel(apertureInfo.apertureVector),2*(apertureInfo.totalNumOfShapes-1));
    
    jacobStruct_dao = padarray(jacobStruct_dao,[0 apertureInfo.totalNumOfShapes-1],0,'post');
    jacobStruct_dos = padarray(jacobStruct_dos,[0 apertureInfo.totalNumOfShapes-1],0,'post');
    % concatenate
    jacobStruct = [jacobStruct_dao; jacobStruct_lfspd; jacobStruct_dosrt; jacobStruct_dos];
end


