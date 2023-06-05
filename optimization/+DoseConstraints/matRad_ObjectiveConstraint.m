classdef matRad_ObjectiveConstraint < DoseConstraints.matRad_DoseConstraint
    % matRad_ObjectiveConstraint implements a wrapper function that turns 
    % dose objectives into constraints
    %
    % use log sum exp approximation, see appendix A in
    % http://scitation.aip.org/content/aapm/journal/medphys/41/8/10.1118/1.4883837
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2020 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Constant)
        name = 'Objective Constraint';
        parameterTypes = {'objFunc'};
        parameterNames = {'f^{max}'};    
    end

    properties
        objective;
        parameters = {1e-5};
        epsilon = 1e-3; %
    end
    
    methods (Access = public)
        function constr = matRad_ObjectiveConstraint(objective,maxObj,epsilon)
            
            %check if objective is a struct and a DoseObjective or Constraint (for init from constraint)
            if isstruct(objective) && ~isempty(strfind(objective.className,'DoseObjectives'))
                objective =  matRad_DoseOptimizationFunction.createInstanceFromStruct(objective);
            end

            if nargin == 1 && isstruct(objective)
                initFromStruct = true;
                inputStruct = objective;
            else
                initFromStruct = false;
                inputStruct = [];
            end

            constr@DoseConstraints.matRad_DoseConstraint(inputStruct);
            
            
            if ~initFromStruct
                
                if nargin == 3 && isscalar(epsilon)
                    constr.epsilon = epsilon;
                end

                if nargin >= 2 && isscalar(maxObj)
                    constr.parameters{1} = maxObj;
                end
                
                if nargin >= 1
                    constr.objective = objective;
                end

            end
            %}
        end
        
        function s = struct(constr)
            s = struct@DoseConstraints.matRad_DoseConstraint(constr);
            s.epsilon = constr.epsilon;
            s.objective = constr.objective;
        end
        
        function cu = upperBounds(constr,n)
            cu = constr.parameters{1}+constr.epsilon;
            %cu = [Inf; constr.parameters{2}];
        end
        function cl = lowerBounds(constr,n)          
            cl = 0;
        end
                
        %% Calculates the Constraint Function value
        function cDose = computeDoseConstraintFunction(constr,dose)
            cDose = constr.objective.computeDoseObjectiveFunction(dose);
        end
        
        %% Calculates the Constraint jacobian
        function cDoseJacob  = computeDoseConstraintJacobian(constr,dose)
            cDoseJacob = constr.objective.computeDoseObjectiveGradient(dose);
        end
        
        function doseParams = getDoseParameters(constr)
            % get only the dose related parameters.
            ix = cellfun(@(c) isequal('dose',c),constr.objective.parameterTypes);
            doseParams = [constr.objective.parameters{ix}];
        end
        
        function constr = setDoseParameters(constr,doseParams)
            % set only the dose related parameters.
            ix = cellfun(@(c) isequal('dose',c),constr.objective.parameterTypes);
            constr.objective.parameters(ix) = num2cell(doseParams);
         
        end        
    

    end
    
end


