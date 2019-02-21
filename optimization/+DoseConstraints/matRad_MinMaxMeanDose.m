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
        function obj = matRad_MinMaxMeanDose(minMeanDose,maxMeanDose)
            
            %If we have a struct in first argument
            if nargin == 1 && isstruct(minMeanDose)
                initFromStruct = true;
                inputStruct = minMeanDose;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            % Call Superclass Constructor (for struct initialization)
            obj@DoseConstraints.matRad_DoseConstraint(inputStruct);
            
            % now handle initialization from other parameters
            if ~initFromStruct
                
                if nargin == 2 && isscalar(maxMeanDose)
                    obj.parameters{2} = maxMeanDose;
                end
                
                if nargin >= 1 && isscalar(minMeanDose)
                    obj.parameters{1} = minMeanDose;
                end
                
            end
        end
        
        function cu = upperBounds(obj,n)
            cu = obj.parameters{2};
        end
        function cl = lowerBounds(obj,n)
            cl = obj.parameters{1};
        end
        
        %Overloads the struct function to add constraint specific
        %parameters
        function s = struct(obj)
            s = struct@DoseConstraints.matRad_DoseConstraint(obj);
            %Nothing to do here...
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


