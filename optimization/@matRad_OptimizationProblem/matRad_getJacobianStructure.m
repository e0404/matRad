function jacobStruct = matRad_getJacobianStructure(optiProb,w,dij,cst)	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% matRad IPOPT callback: jacobian structure function for inverse planning supporting max dose	
% constraint, min dose constraint, min mean dose constraint, max mean dose constraint, 	
% min EUD constraint, max EUD constraint, max DVH constraint, min DVH constraint 	
% 	
% call	
%   jacobStruct = matRad_getJacobStruct(dij,cst)	
%	
% input	
%   dij: dose influence matrix	
%   cst: matRad cst struct	
%	
% output	
%   jacobStruct: jacobian of constraint function	
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
jacobStruct = sparse([]);	
 % compute objective function for every VOI.	
for i = 1:size(cst,1)	
     % Only take OAR or target VOI.	
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )	
         % loop over the number of constraints for the current VOI	
        for j = 1:numel(cst{i,6})	
            	
            obj = cst{i,6}{j};	
            	
            % only perform computations for constraints	
              if isa(obj,'DoseConstraints.matRad_DoseConstraint')	
                	
                % get the jacobian structure depending on dose	
                jacobDoseStruct = obj.getDoseConstraintJacobianStructure(numel(cst{i,4}{1}));	
                nRows = size(jacobDoseStruct,2);	
                jacobStruct = [jacobStruct; repmat(spones(mean(dij.physicalDose{1}(cst{i,4}{1},:))),nRows,1)];	
                 
             end	
         end	
     end	
 end