function c = matRad_constraintFunctions(optiProb,I,dij,w,cst)
% matRad IPOPT callback: constraint function for inverse planning 
% supporting max dose constraint, min dose constraint, min mean dose constraint, 
% max mean dose constraint, min EUD constraint, max EUD constraint, 
% max DVH constraint, min DVH constraint 
% 
% call
%   c = matRad_constraintFunctions(optiProb,w,dij,cst)
%
% input
%   optiProb:   option struct defining the type of optimization
%   I:          intensity vector
%   dij:        dose influence matrix
%   cst:        matRad cst struct
%
% output
%   c:          value of constraints
%
% References
%   -
%
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


% get current dose / effect / RBExDose vector
%d = matRad_backProjection(w,dij,optiProb);
optiProb.BP = optiProb.BP.compute(dij,w);
dCombined = optiProb.BP.GetResult();
index = dij.doseGrid.numOfVoxels;
d{1} = dCombined{1}(1:index);
DADR{1} = dCombined{1}(index+1:end);

% Initializes constraints
c = [];

% compute objective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR 
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') ||  ( isequal(cst{i,3},'TARGET') ))

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
            obj = cst{i,6}{j};
            
            if isa(obj,'DoseConstraints.matRad_DoseConstraint')      
                d_i = d{1}(cst{i,4}{1});
                c = [c; obj.computeDoseConstraintFunction(d_i)];
            elseif isa(obj,'DADRConstraints.matRad_DADRConstraint')  
                dadr_i = DADR{1}(cst{i,4}{1});
                c = [c; obj.computeDoseConstraintFunction(dadr_i)];
            end
            
                %end
            
        end

    end % if structure not empty and oar or target

end % over all structures