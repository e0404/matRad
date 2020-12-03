classdef matRad_OptimizationProblem < handle
    %handle class to keep state easily
    
    properties
        BP
        bioOpt = '';
        
        auxVars = [];
    end
    
    %properties (Access = private)
    %    currentDose
    %    currentWeights
    %end
    
    methods
        function obj = matRad_OptimizationProblem(backProjection)
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
        
        %lagrangian hessian declaration
        lHessian = matRad_lagrangianHessian(optiProb,w,dij,cst,lambda,sigma);
        
        %Jacobian Structure
        jacobStruct = matRad_getJacobianStructure(optiProb,w,dij,cst)
        
        %Hessian Structure
        hessianStruct = matRad_getHessianStructure(optiProb,w,dij,cst);              
        
        [cl,cu] = matRad_getConstraintBounds(optiProb,cst);
        
        %Auxiliary Variables
        function initializeAuxiliaryVariables(optiProb,cst)
            z = exp((dose-max(dose))/obj.epsilon);
            
        end
        
        function lb = lowerBounds(optiProb,w)
            lb = zeros(size(w));
        end
        
        function ub = upperBounds(optiProb,w)
            ub = Inf * ones(size(w));
        end
    end
end

