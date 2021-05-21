classdef matRad_OptimizationProblem < handle
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
        BP
        quantityOpt = '';
        useMaxApprox = 'logsumexp'; %'pnorm'; %'logsumexp'; %'none';
        p = 30; %Can be chosen larger (closer to maximum) or smaller (closer to mean). Only tested 20 >= p >= 1
    end
    
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
    
    methods (Access = protected)
        function [val,grad] = logSumExp(optiProb,fVals)
            % https://doi.org/10.1093/imanum/draa038
            [fMax,ixMax] = max(fVals(:));
            
            ix = true(numel(fVals),1);
            ix(ixMax) = 0;

            tmp = exp(fVals - fMax);
                       
            expSum = sum(tmp(ix));
            val = fMax + log1p(expSum); %log1p(x) = Matlab's numerically accurate log(1+x) 
            
            grad = tmp ./ (1 + expSum);
        end
        
        function [val,grad] = pNorm(optiProb,fVals,n)
            if nargout < 3
                n = numel(fVals);
            end
            
            p = optiProb.p;
            
            valMax = max(fVals(:));            
            tmp = fVals./valMax;            
            
            pNormVal = sum(tmp(:).^p)^(1/p);
            
            fac = (1/n)^(1/p);
            val = valMax*fac*pNormVal;

            grad = fac * (tmp ./ pNormVal).^(p-1);
        end
    end                
end

