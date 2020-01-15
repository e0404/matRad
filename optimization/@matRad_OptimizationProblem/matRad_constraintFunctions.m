function c = matRad_constraintFunctions(optiProb,w,dij,cst)

% matRad IPOPT callback: constraint function for inverse planning supporting max dose
% constraint, min dose constraint, min mean dose constraint, max mean dose constraint,
% min EUD constraint, max EUD constraint, max DVH constraint, min DVH constraint 
% 
% call
%   jacob = matRad_jacobFunc(w,dij,cst,optiProb)
%
% input
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   cst:     matRad cst struct
%   optiProb: option struct defining the type of optimization
%
% output
%   c: value of constraints
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
d = optiProb.BP.GetResult();

% Initializes constraints
c = [];

% compute objective function for every VOI.
for  i = 1:size(cst,1)

    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
            
            obj = cst{i,6}{j};
            
            % only perform computations for constraints
            % if ~isempty(strfind(obj.type,'constraint'))
            if isa(obj,'DoseConstraints.matRad_DoseConstraint')
                
                % if we have effect optimization, temporarily replace doses with effect
                % Maybe we should put some switch into the classes for that
                if (~isequal(obj.name, 'max dose constraint')      && ~isequal(obj.name, 'min dose constraint')      &&...
                    ~isequal(obj.name, 'max mean dose constraint') && ~isequal(obj.name, 'min mean dose constraint') && ...
                    ~isequal(obj.name, 'min EUD constraint')       && ~isequal(obj.name, 'max EUD constraint'))     && ...
                    (isa(optiProb.BP,'matRad_EffectProjection') && ~isa(optiProb.BP,'matRad_VariableRBEProjection'))
                    
                    doses = obj.getDoseParameters();
                    
                    effect = cst{i,5}.alphaX*doses + cst{i,5}.betaX*doses^2;
                    
                    obj = obj.setDoseParameters(effect);
                end

                % if conventional opt: just add constraints of nominal dose
                %if strcmp(cst{i,6}(j).robustness,'none')

                    d_i = d{1}(cst{i,4}{1});

                    %c = [c; matRad_constFunc(d_i,cst{i,6}(j),d_ref)];
                    c = [c; obj.computeDoseConstraintFunction(d_i)];
                    
                %else
                    
                    %error('Invalid robustness setting.');

                %end % if we are in the nominal sceario or rob opt
            
            end % if we are a constraint

        end % over all defined constraints & objectives

    end % if structure not empty and oar or target

end % over all structures
