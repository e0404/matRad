classdef matRad_MinMaxMeanDose < DoseConstraints.matRad_DoseConstraint
    %MATRAD_DOSEOBJECTIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        name = 'mean dose constraint';
        parameterNames = {'\mu_d^{min}', '\mu_d^{max}'};
        %parameterIsDose = logical([1 1]);
        parameterTypes = {'dose','dose'};
    end
    
    properties
        parameters = {0,30};
    end
        
    methods
        function cu = upperBounds(obj,n)
            cu = obj.parameters{2};
        end
        function cl = lowerBounds(obj,n)
            cl = obj.parameters{1};
        end
        %% Calculates the Constraint Function value
        function cDose = computeDoseConstraintFunction(obj,dose)                       
            cDose = mean(dose);        
        end
        
        %% Calculates the Constraint jacobian
        function cDoseJacob  = computeDoseConstraintJacobian(obj,dose)
            cDoseJacob = ones(numel(dose),1)./numel(dose);
        end
    end
    
end


