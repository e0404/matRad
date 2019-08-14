classdef matRad_OptimizationProblem < handle
    %handle class to keep state easily
    
    properties
        BP
        bioOpt = '';
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
        
        %Jacobian Structure
        jacobStruct = matRad_getJacobianStructure(optiProb,w,dij,cst)
        
        [cl,cu] = matRad_getConstraintBounds(optiProb,cst)
        
        function lb = lowerBounds(optiProb,w)
            lb = zeros(size(w));
        end
        
        function ub = upperBounds(optiProb,w)
            ub = Inf * ones(size(w));
        end
    end
end

