classdef matRad_OptimizationProblem < handle
    %handle class to keep state easily
    
    properties
        BP %backProjection
        bioOpt = '';
    end
    
    properties (Access = private)
        currentDose
        currentWeights
    end
    
    methods
        function obj = matRad_OptimizationProblem(backProjection,dij,cst)
            obj.BP = backProjection;
        end       
        
        %Objective function declaration
        fVal = matRad_objectiveFunction(optiProb,w,dij,cst)   
        
        %Objective gradient declaration
        fGrad = matRad_objectiveGradient(optiProb,w,dij,cst)
        
        %Constraint function declaration
        cVal = matRad_constraintFunctions(optiProb,w,dij,cst)
        
        %Constraint Jacobian declaration
        cJacob = matRad_constraintJacobian(optiProb,w,dij,cst)
        
        [cl,cu] = matRad_getConstraintBounds(optiProb,cst)
        
        function lb = lowerBounds(optiProb)
            lb = zeros(size(optiProb.currentWeights));
        end
        
        function ub = upperBounds(optiProb)
            ub = Inf * ones(size(optiProb.currentWeights));
        end
    end
end

