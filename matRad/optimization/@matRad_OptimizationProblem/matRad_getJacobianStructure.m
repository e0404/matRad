function jacobStruct = matRad_getJacobianStructure(optiProb,w,dij,cst)	
% matRad IPOPT callback: jacobian structure function for inverse planning 
% supporting max dose constraint, min dose constraint, min mean dose constraint, 
% max mean dose constraint, min EUD constraint, max EUD constraint, max DVH 
% constraint, min DVH constraint 	
% 	
% call	
%   jacobStruct = matRad_getJacobStruct(optiProb,w,dij,cst)	
%	
% input	
%   optiProb: matRad optimization problem
%   w:        beamlet/ pencil beam weight vector
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
for i = 1:size(optiProb.constrIdx,1)	
   obj = optiProb.constraints{i};
   curConIdx = optiProb.constrIdx(i,1);

   % get the jacobian structure depending on dose	
   jacobDoseStruct = obj.getDoseConstraintJacobianStructure(numel(cst{curConIdx,4}{1}));	
   nRows = size(jacobDoseStruct,2);	
   jacobStruct = [jacobStruct; repmat(spones(mean(dij.physicalDose{1}(cst{curConIdx,4}{1},:))),nRows,1)];	
      
end