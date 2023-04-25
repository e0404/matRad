classdef matRad_OptimizationProblemMixMod < handle
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
        objectives = {}; %cell array storing all objectives, has to be initialized at the start
        constraints = {}; %
        objidx;
        constridx;
        normalizationScheme = struct('type','none');

        useMaxApprox = 'logsumexp'; %'pnorm'; %'logsumexp'; %'none';
        p = 30; %Can be chosen larger (closer to maximum) or smaller (closer to mean). Only tested 20 >= p >= 1

        minimumW = NaN;
        maximumW = NaN;
    end
    
    methods
        function obj = matRad_OptimizationProblemMixMod(backProjection,cst)
        
            obj.BP = backProjection; %needs to be initalized to have access to setBiologicalDosePrescriptions
                if nargin == 2
                    obj.extractObjectivesAndConstraintsFromcst(cst);
                end
        end       
        
        %Objective function declaration
        fVal = matRad_objectiveFunction2(optiProb,w,dij,cst)   
        %fVal = matRad_objectiveFunction(optiProb,w,dij,cst) 
        
        %Objective gradient declaration
        fGrad = matRad_objectiveGradient(optiProb,w,dij,cst)
        %fGrad = matRad_objectiveGradient2(optiProb,w,dij,cst)

        %Constraint function declaration
        cVal = matRad_constraintFunctions(optiProb,w,dij,cst)
        
        %Constraint Jacobian declaration
        cJacob = matRad_constraintJacobian(optiProb,w,dij,cst)
        
        %Jacobian Structure
        jacobStruct = matRad_getJacobianStructure(optiProb,w,dij,cst)
        
        [cl,cu] = matRad_getConstraintBounds(optiProb,cst)
        
        function lb = lowerBounds(optiProb,w)
            minW = optiProb.minimumW;
            if isnan(minW)
                lb = zeros(size(w));
            elseif isscalar(minW)
                lb = minW*ones(size(w));
            elseif isvector(minW) && all(size(minW) == size(w))
                lb = minW;
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Minimum Bounds for Optimization Problem could not be set!');
            end
        end
        
        function ub = upperBounds(optiProb,w)
            maxW = optiProb.maximumW;
            if isnan(maxW)
                ub = Inf(size(w));
            elseif isscalar(maxW)
                ub = maxW*ones(size(w));
            elseif isvector(maxW) && all(size(maxW) == size(w))
                ub = maxW;
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Maximum Bounds for Optimization Problem could not be set!');
            end
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
    end          

    methods (Access = private)
        function extractObjectivesAndConstraintsFromcst(optiProb,cst)
            %used to store objectives in cell array as property of optimization Problem
            optiProb.objidx = [];
            optiProb.constridx = [];
            optiProb.objectives = {};
            optiProb.constraints = {};
            
            for i = 1:size(cst,1) % loop over cst
                
                if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

                    for j = 1:numel(cst{i,6})
                        %check whether dose objective or constraint
                        obj = cst{i,6}{j};
                        if isstruct(cst{i,6}{j})
                            obj =  matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
                        end
                        if contains(class(obj),'DoseObjectives')
                            optiProb.objidx = [optiProb.objidx;i,j];
                            obj = optiProb.BP.setBiologicalDosePrescriptions(obj,cst{i,5}.alphaX,cst{i,5}.betaX);
                            optiProb.objectives(end+1) = {obj};

                        elseif contains(class(obj),'DoseConstraints')
                            optiProb.constridx = [optiProb.constridx;i,j];
                            obj = optiProb.BP.setBiologicalDosePrescriptions(obj,cst{i,5}.alphaX,cst{i,5}.betaX);
                            optiProb.constraints(end+1) = {obj};
                        end
                    end
                end
            end
        end


    end        
end

