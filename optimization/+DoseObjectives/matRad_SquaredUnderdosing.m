classdef matRad_SquaredUnderdosing < DoseObjectives.matRad_DoseObjective
    %MATRAD_DOSEOBJECTIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        name = 'Squared Underdosing';
        parameterNames = {'d^{min}'};
        parameterTypes = {'dose'};
    end
    
    properties
        parameters = {60};
        penalty = 1;
    end
    
    methods 
        %% Calculates the Objective Function value
        function fDose = computeDoseObjectiveFunction(obj,dose)                       
            % overdose : dose minus prefered dose
            underdose = dose - obj.parameters{1};
            
            % apply positive operator
            underdose(underdose>0) = 0;

            % claculate objective function
            fDose = obj.penalty/numel(dose) * (underdose'*underdose);
        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad   = computeDoseObjectiveGradient(obj,dose)
            % overdose : dose minus prefered dose
            underdose = dose - obj.parameters{1};
            
            % apply positive operator
            underdose(underdose>0) = 0;

            % calculate delta
            fDoseGrad = 2 * obj.penalty/numel(dose) * underdose;
        end
    end
    
end
