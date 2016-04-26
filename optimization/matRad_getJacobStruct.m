function jacobStruct = matRad_getJacobStruct(dij,cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: jacobian structure function for inverse planning supporting max dose
% constraint, min dose constraint, min max dose constraint, min mean, max
% min, min max mean constraint, min EUD constraint, max EUDconstraint, 
% min max EUD constraint, exact DVH constraint, max DVH constraint, 
% min DVH constraint 
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

% Initializes constraints
jacobStruct = sparse([]);

% compute objective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
                    
        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
                        
            if isequal(cst{i,6}(j).type, 'max dose constraint') || ...
               isequal(cst{i,6}(j).type, 'min dose constraint') || ...
               isequal(cst{i,6}(j).type, 'max mean dose constraint') || ...
               isequal(cst{i,6}(j).type, 'min mean dose constraint') || ...
               isequal(cst{i,6}(j).type, 'min max mean dose constraint') || ...
               isequal(cst{i,6}(j).type, 'max EUD constraint') || ...
               isequal(cst{i,6}(j).type, 'min EUD constraint') || ...
               isequal(cst{i,6}(j).type, 'min max EUD constraint') || ...
               isequal(cst{i,6}(j).type, 'exact DVH constraint') || ...
               isequal(cst{i,6}(j).type, 'max DVH constraint') || ... 
               isequal(cst{i,6}(j).type, 'min DVH constraint')

                jacobStruct = [jacobStruct; spones(mean(dij.physicalDose(cst{i,4},:)))];
                               
            end
      
        end
        
    end
end
   
end