classdef matRad_DoseConstraintFromObjective < DoseConstraints.matRad_DoseConstraint
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
        parameterTypes = {'objFunc','scalar'};
        parameterNames = {'f^{max}','slackParameter'};    
    end

    properties
        objective;
        parameters = {1e-5, 1e-3};
    end
    
    methods (Access = public)
        function this = matRad_DoseConstraintFromObjective(objective,maxObj,slackParameter)
            
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

            this@DoseConstraints.matRad_DoseConstraint(inputStruct);
            
            
            if ~initFromStruct
                
                if nargin == 3 && isscalar(slackParameter)
                    this.parameters{2} = slackParameter;
                end

                if nargin >= 2 && isscalar(maxObj)
                    this.parameters{1} = maxObj;
                end
                
                if nargin >= 1
                    this.objective = objective;
                end

            end
            %}
        end
        
        function s = struct(this)
            s = struct@DoseConstraints.matRad_DoseConstraint(this);
            s.objective = this.objective;
        end
        
        function cu = upperBounds(this,n)
            cu = this.parameters{1}+this.slackParameter;
            %cu = [Inf; this.parameters{2}];
        end
        function cl = lowerBounds(this,n)          
            cl = 0;
        end
                
        %% Calculates the Constraint Function value
        function cDose = computeDoseConstraintFunction(this,dose)
            cDose = this.objective.computeDoseObjectiveFunction(dose);
        end
        
        %% Calculates the Constraint jacobian
        function cDoseJacob  = computeDoseConstraintJacobian(this,dose)
            cDoseJacob = this.objective.computeDoseObjectiveGradient(dose);
        end
        
        function doseParams = getDoseParameters(this)
            % get only the dose related parameters.
            ix = cellfun(@(c) isequal('dose',c),this.objective.parameterTypes);
            doseParams = [this.objective.parameters{ix}];
        end
        
        function this = setDoseParameters(this,doseParams)
            % set only the dose related parameters.
            ix = cellfun(@(c) isequal('dose',c),this.objective.parameterTypes);
            this.objective.parameters(ix) = num2cell(doseParams);
         
        end        
    

    end
    
end


