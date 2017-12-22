classdef matRad_SquaredDeviation < DoseObjectives.matRad_DoseObjective
    %MATRAD_DOSEOBJECTIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        name = 'Squared Deviation'
    end
    
    properties
        parameters = {'d^{ref}'; 60}
        penalty = 1
    end
    
    methods 
        %% Calculates the Objective Function value
        function fDose = computeDoseObjectiveFunction(obj,dose)
            % deviation : dose minus prefered dose
            deviation = dose - obj.parameters{2,1};
            % claculate objective function
            fDose = obj.penalty/numel(dose) * (deviation'*deviation);
        end
        
        %% Calculates the Objective Function Gradient
        function fDoseGrad   = computeDoseObjectiveGradient(obj,dose)
            % deviation : Dose minus prefered dose
            deviation = dose - obj.parameters{2,1};

            % calculate delta
            fDoseGrad = 2 * obj.penalty/numel(dose) * deviation;
        end
    end
    
end
