classdef matRad_OptimizationProblemDADR < matRad_OptimizationProblem
    %matRad_OptimizationProblem Main Class for fluence optimization problems
    % Describes a standard fluence optimization problem by providing the 
    % implementation of the objective & constraint function/gradient wrappers
    % and managing the mapping and backprojection of the respective dose-
    % related quantity
    %
    % References
    %   -
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2020 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        BP_dadr              %matRad_BackProjection object for mapping & backprojection        
    end
    
    methods
        function obj = matRad_OptimizationProblemDADR(doseProjection,dadrProjection)
            obj@matRad_OptimizationProblem(doseProjection);
            obj.BP_dadr = dadrProjection;
        end       
        
        %Objective function declaration
        [fVal] = matRad_objectiveFunction(optiProb,w,dij,cst)   
        
        %Objective gradient declaration
        [fGrad] = matRad_objectiveGradient(optiProb,w,dij,cst)
        
        %Constraint function declaration
        cVal = matRad_constraintFunctions(optiProb,w,dij,cst)
        
        %Constraint Jacobian declaration
        cJacob = matRad_constraintJacobian(optiProb,w,dij,cst)
        
        %Jacobian Structure
        jacobStruct = matRad_getJacobianStructure(optiProb,I,dij,cst)
        
        [cl,cu] = matRad_getConstraintBounds(optiProb,cst)
       
        
        function lb = lowerBounds(~,wCombined)
            %lb = 2.9*10^3*ones(size(I));  
            lb = zeros(size(wCombined));
        end
        
        function ub = upperBounds(~,wCombined)
            %             % Include upper bounds for the intensity
            %ub = 2.9*10^6*ones(size(I)); %arbitrary value for now
            ub = Inf*ones(size(wCombined));
        end
    end
end

