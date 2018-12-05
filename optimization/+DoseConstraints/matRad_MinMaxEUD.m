classdef matRad_MinMaxEUD < DoseConstraints.matRad_DoseConstraint
    %MATRAD_DOSEOBJECTIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        name = 'EUD constraint';
        parameterNames = {'k','EUD^{min}', 'EUD^{max}'};
        %parameterIsDose = logical([0 1 1]);
        parameterTypes = {'numeric','dose','dose'};
    end
    
    properties
        parameters = {5,0,30};
    end
        
    methods
        function cu = upperBounds(obj,n)
            cu = obj.parameters{3};
        end
        function cl = lowerBounds(obj,n)
            cl = obj.parameters{2};
        end
        %% Calculates the Constraint Function value
        function cDose = computeDoseConstraintFunction(obj,dose)                       
            k = obj.parameters{1};
            cDose = mean(dose.^k)^(1/k);           
        end
        
        %% Calculates the Constraint jacobian
        function cDoseJacob  = computeDoseConstraintJacobian(obj,dose)
            k = obj.parameters{1};
            cDoseJacob = nthroot(1/numel(dose),k) * sum(dose.^k)^((1-k)/k) * (dose.^(k-1));
        end
    end
    
end


