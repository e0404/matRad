classdef matRad_DoseMiniMax < DoseObjectives.matRad_DoseObjective
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    %
    properties (Constant)
        name = 'Dose MiniMax Objective';         	
        parameterNames	= {'d^{max}'};
        parameterTypes	= {'dose'};       
    end
    properties
        penalty = 1;
        parameters = {30};
        epsilon = 1e-3; %slack parameter for the logistic approximation
    end
    
    methods
        function obj = matRad_DoseMiniMax(penalty,dRef)
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
                if nargin == 2 && isscalar(dRef)
                    obj.parameters{1} = dRef;
                end

                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
        
        function s = struct(obj)
            s = struct@matRad_DoseOptimizationFunction(obj);
        end
        
        function f = computeDoseObjectiveFunction(obj,dose)
            %auxVar holds the current maximum dose
            %auxVar * obj.penalty - obj.parameters{1};
            %f = max([max(dose) - obj.parameters{1} , 0]) * obj.penalty;
            f = max([max(dose) - obj.parameters{1}, 0]) * obj.penalty;
        end
        
        function fGrad = computeDoseObjectiveGradient(obj,dose)
            fGrad = dose - max(dose);
            ix = dose < obj.parameters{1};
            fGrad = exp(fGrad / obj.epsilon);
            fGrad = obj.penalty/sum(fGrad) * fGrad;
            fGrad(ix) = 0;
        end
        
        function fHess = computeDoseObjectiveHessian(obj,dose)
            
        end
        
        %{
        %% Implement the constraint part
        function cu = upperBounds(obj,n)
            cu = Inf;
        end
        function cl = lowerBounds(obj,n)
            cl = 0;
        end
        %}
        
        %{
        %Calcualte the constraint function
        function cDose = computeDoseConstraintFunction(obj,dose)
            %dMax = max(dose);
            %cDose = dose - auxVar; 
            cDose = dose - max(dose);
        end
        %}
        
        
        
        %{
        %% Calculates the Constraint jacobian
        function cDoseJacob  = computeDoseConstraintJacobian(obj,dose)
            cDoseJacob = speye(numel(dose),numel(dose));
        end
        
        function cDoseHessian = computeDoseConstraintHessian(obj,dose,lambda)
            cDoseHessian = sparse(numel(dose),numel(dose));
        end
        %}
        
    end
end

