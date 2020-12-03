function hessianStruct = matRad_getHessianStructure(optiProb,w,dij,cst)	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% matRad IPOPT callback: hessian structure function for inverse planning supporting max dose	
% constraint, min dose constraint, min mean dose constraint, max mean dose constraint, 	
% min EUD constraint, max EUD constraint, max DVH constraint, min DVH constraint 	
% 	
% call	
%   hessianStruct = matRad_getHessianStructure(dij,cst)	
%	
% input	
%   dij: dose influence matrix	
%   cst: matRad cst struct	
%	
% output	
%   hessianStruct: sparse structure of hessian of the lagrangian 	
%	
% References	
%	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%	
% Copyright 2016 the matRad development team. 	
% 	
% This file is part of the matRad project. It is subject to the license 	
% terms in the LICENSE file found in the top-level directory of this 	
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 	
% of the matRad project, including this file, may be copied, modified, 	
% propagated, or distributed except according to the terms contained in the 	
% LICENSE file.	
%	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
 % Initializes constraints	
%hessianStructDose = sparse(dij.numVoxels,dij.numVoxels);
allIx = [];
 % compute objective function for every VOI.	
for i = 1:size(cst,1)	
     % Only take OAR or target VOI.	
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )	
         % loop over the number of constraints for the current VOI	
        for j = 1:numel(cst{i,6})	
            	
            obj = cst{i,6}{j};	
            voiIx = cst{i,4}{1};
            allIx = [allIx; voiIx];
            	
            %if isa(obj,'DoseConstraints.matRad_DoseConstraint')	
            %    spStruct = sparse(obj.getDoseConstraintHessianStructure(numel(voiIx)));
            %else
            %    spStruct = sparse(obj.getDoseObjectiveHessianStructure(numel(voiIx)));
            %end
            
            %hessianStructDose(voiIx,voiIx) = hessianStructDose(voiIx,voiIx) + spStruct;
            
         end	
     end	
end
 
%All elements with structure to ones
%hessianStructDose(hessianStructDose >= 1) = 1;

allIx = unique(allIx);

hessianStruct = spones(mean(dij.physicalDose{1}(allIx,:)));	
hessianStruct = hessianStruct' * hessianStruct;

hessianStruct = tril(hessianStruct);




