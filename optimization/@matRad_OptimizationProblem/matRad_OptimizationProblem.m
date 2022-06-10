classdef matRad_OptimizationProblem < handle
    %matRad_OptimizationProblem Main Class for fluence optimization problems
    % Describes a standard fluence optimization problem by providing the
    % implementation of the objective & constraint function/gradient wrappers
    % and managing the mapping and backprojection of the respective dose-
    % related quantity
    %
    % References
    %   [1] https://doi.org/10.1093/imanum/draa038
    %   [2] https://doi.org/10.1002/mp.14148
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
        
        % Maximum approximation
        useMaxApprox = 'cheapCOWC'; %'pnorm'; %'cheapCOWC'; %'logsumexp'; %'none';
        
        % Parameters for pnorm approximation
        p = 30; %Can be chosen larger (closer to maximum) or smaller (closer to mean). Only tested 20 >= p >= 1
        
        % Paremeters for CheapCOWC
        p1=1; %
        p2=4; %
        
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
            % [1] stable log sum exp trick
            [fMax,ixMax] = max(fVals(:));
            
            ix = true(numel(fVals),1);
            ix(ixMax) = 0;
            
            tmp = exp(fVals - fMax);
            
            expSum = sum(tmp(ix));
            val = fMax + log1p(expSum); %log1p(x) = Matlab's numerically accurate log(1+x)
            
            grad = tmp ./ (1 + expSum);
        end
        
        function [val,grad] = pNorm(optiProb,fVals,n)
            % Implemented as proposed in [2] including a normalization for stability of the exponent.
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
        
        % Cheap-minimax gradient calculation 
        function [val,grad] = cheapCOWC(optiProb,fVals,fProb)
            
            mask = true(size(fVals));
            val_array = zeros(numel(fVals),1);
            grad_array = cell(numel(fVals),1);
            norm_array = zeros(numel(fVals),1);
            ix = true(numel(fVals),1);
            
            
            for i = 1: numel(fVals)
                [fMax,ixMax] = max(fVals(:));
                
                ix(ixMax) = 0;
                
                tmp = exp(fVals - fMax);
                
                expSum = sum(tmp(ix));
                
                val_array(i) = (fMax + log1p(expSum))*fProb(ixMax);
                grad_array{i} = tmp.*fProb(ixMax).*mask;
                norm_array(i) = (1 + expSum).*fProb(ixMax);
                
                fVals(ixMax) = 0;
                fProb(ixMax) = 0;
                mask(ixMax) = false;
            end
            
            val = 0;
            grad = zeros(size(fVals));
            norm = sum(norm_array(optiProb.p1:optiProb.p2));
            
            for i = optiProb.p1:optiProb.p2
                val = val + val_array(i);
                grad = grad + grad_array{i}/norm;
            end
            
        end
        
    end
end

