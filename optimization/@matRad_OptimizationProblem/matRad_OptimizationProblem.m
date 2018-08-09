classdef matRad_OptimizationProblem < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        BP
        bioOpt = '';
    end
    
    properties (Access = private)
        currentDose
        currentWeights
    end
    
    methods
        function obj = matRad_OptimizationProblem(matRad_BackProjection,dij,cst)
            obj.BP = matRad_BackProjection;
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

