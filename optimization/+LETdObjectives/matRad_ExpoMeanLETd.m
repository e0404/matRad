classdef matRad_ExpoMeanLETd < LETdObjectives.matRad_LETdObjective
% matRad_ExpoMeanLETd implements a penalized exponential mean LETd objective
% f = ((sum(LETd)/N)-LETdRef)^k
%   See matRad_LETdObjective for interface description
%
% References
%   matRad_MeanDose
%   matRad_SquaredDeviation
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
    
    properties (Constant)
        name = 'Exponential Mean LETd';
        parameterNames = {'LETd^{ref}','k'};       % k is the exponent of the mean function
        parameterTypes = {'LETd','numeric'};
    end
    
    properties
        parameters = {0,2};        
        penalty = 1;
    end
    
    methods 
        function obj = matRad_ExpoMeanLETd(penalty,LETdMeanRef,k)
           
            % if we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            % Call Superclass Constructor (for struct initialization)
            obj@LETdObjectives.matRad_LETdObjective(inputStruct);
            
            % now handle initialization from other parameters
            if ~initFromStruct
                if nargin == 3 && isscalar(k)
                    obj.parameters{2} = k;
                end
                if nargin >= 2 && isscalar(LETdMeanRef)
                    obj.parameters{1} = LETdMeanRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
                   
            %% Downwards compatability / set default values
            % TODO: maybe move into set method for parameters
            if numel(obj.parameters) < 1
                obj.parameters{1} = 0;
            end

            if numel(obj.parameters) < 2
                obj.parameters{2} = 1;
            end
        end
        %% Calculates the Objective Function value
        function fLETd = computeLETdObjectiveFunction(obj,LETd)
            fLETd =  obj.objectiveExpoDiff(LETd);
        end
        
        %% Calculates the Objective Function gradient
        function fLETdGrad   = computeLETdObjectiveGradient(obj,LETd)
            fLETdGrad = obj.gradientExpoDiff(LETd);
        end
    end
     methods (Access = protected)
        function fLETd = objectiveExpoDiff(obj,LETd)
            % Formula for the Objective Function value
            fLETd = (mean(LETd(:)) - obj.parameters{1})^obj.parameters{2}; 
        end

        function fLETdGrad = gradientExpoDiff(obj,LETd)
            % Derivative of the Objective Function
            fLETdGrad = obj.parameters{2}*((mean(LETd(:))-obj.parameters{1}))^(obj.parameters{2}-1) * ones(size(LETd(:)))/numel(LETd);
        end

    end
    
end


