classdef matRad_MeanDose < DoseObjectives.matRad_DoseObjective
    %MATRAD_MEANDOSE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        name = 'Mean Dose';
        parameterNames = {'d^{ref}',{'difference'}};
        parameterTypes = {'dose',{'linear','squared'}};
    end
    
    properties
        parameters = {0,1};        
        penalty = 1;
    end
    
    methods 
        function obj = matRad_MeanDose(penalty,dMeanRef)
           
            % if we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DoseObjectives.matRad_DoseObjective(inputStruct);
            
            if ~initFromStruct
                if nargin == 2 && isscalar(dMeanRef)
                    obj.parameters{1} = dMeanRef;
                end

                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
            
        end       
        
        %% Calculates the Objective Function value
        function fDose = computeDoseObjectiveFunction(obj,dose)
            %fDose = obj.penalty * abs(mean(dose(:)) - obj.parameters{1});
            if obj.parameters{2} == 2
                fDose = obj.penalty * (mean(dose(:)) - obj.parameters{1})^2;
            else
            	fDose = obj.penalty * (mean(dose(:)) - obj.parameters{1});
            end
        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad   = computeDoseObjectiveGradient(obj,dose)
            
            n = numel(dose);
            
            %fDoseGrad = (obj.penalty/numel(dose))*sign(dose(:)-obj.parameters{1});
            if obj.parameters{2} == 2
                fDoseGrad = obj.penalty*2*(mean(dose(:))-obj.parameters{1}) * ones(n,1)/n;
            else
                fDoseGrad = obj.penalty * ones(size(dose(:))) / numel(dose);
            end
        end
        
        %% Calculates the Objective Function Hessian
        function fDoseHess   = computeDoseObjectiveHessian(obj,dose)         
            
            n = numel(dose);
            
            if obj.parameters{2} == 2
                fDoseHess = 2 * obj.penalty/numel(dose)^2 * ones(n);
            else
                fDoseHess = sparse(numel(dose),numel(dose));
            end
        end
    end
    
end

