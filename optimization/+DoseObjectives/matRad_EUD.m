classdef matRad_EUD < DoseObjectives.matRad_DoseObjective
    %MATRAD_DOSEOBJECTIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        name = 'EUD';
        parameterNames = {'EUD^{ref}', 'k'};
        parameterTypes = {'dose','numeric'};
    end
    
    properties
        parameters = {0, 3.5};
        penalty = 1;
    end
    
    methods
        function obj = matRad_EUD(penalty,eudRef, eudExponent)
            %If we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DoseObjectives.matRad_DoseObjective(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                if nargin >= 3 && isscalar(eudExponent)
                    obj.parameters{2} = eudExponent;
                end
                
                if nargin >= 2 && isscalar(eudRef)
                    obj.parameters{1} = eudRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
        
        %% Calculates the Objective Function value
        function fDose = computeDoseObjectiveFunction(obj,dose)
            % get exponent for EUD
            k = obj.parameters{2};
            
            % calculate power sum
            powersum = sum(dose.^k);
            
            
            
            %Calculate objective
            
            %This check is not needed since dose is always positive
            %if powersum > 0
            fDose = obj.penalty * (nthroot(powersum/numel(dose),k) - obj.parameters{1})^2;
            %end
        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad  = computeDoseObjectiveGradient(obj,dose)
            % get exponent for EUD
            k = obj.parameters{2};
            
            % calculate power sum
            powersum = sum(dose.^k);
            %This check is not needed since dose is always positive
            %if powersum > 0
            
            %derivatives = nthroot(1/numel(dose),k) * powersum^((1-k)/k) * (dose.^(k-1));
            fDoseGrad = 2 * obj.penalty * nthroot(1/numel(dose),k) * powersum^((1-k)/k) * (dose.^(k-1)) .* (nthroot(powersum/numel(dose),k) - obj.parameters{1});
            %end
            if any(~isfinite(fDoseGrad)) % check for inf and nan for numerical stability
                error(['EUD computation failed. Reduce exponent to resolve numerical problems.']);
            end
        end
    end
    
end

