classdef matRad_ExpoMeanLETxDose < LETxDoseObjectives.matRad_LETxDoseObjective
% matRad_ExpoMeanLETxDose implements a penalized exponential mean LETxDose objective
% f = ((sum(LETxDose)/N)-LETxDoseRef)^k
%   See matRad_LETxDoseObjective for interface description
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
        name = 'Exponential Mean LETxDose';
        parameterNames = {'LETxDose^{ref}','k'};       % k is the exponent of the mean function
        parameterTypes = {'LETxDose','numeric'};
    end
    
    properties
        parameters = {0,2};        
        penalty = 1;
    end
    
    methods 
        function obj = matRad_ExpoMeanLETxDose(penalty,LETxDoseRef,k)
           
            % if we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            % Call Superclass Constructor (for struct initialization)
            obj@LETxDoseObjectives.matRad_LETxDoseObjective(inputStruct);
            
            % now handle initialization from other parameters
            if ~initFromStruct
                if nargin == 3 && isscalar(k)
                    obj.parameters{2} = k;
                end
                if nargin >= 2 && isscalar(LETxDoseRef)
                    obj.parameters{1} = LETxDoseRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
                   
            %% Downwards compatability / set default values
            
            if numel(obj.parameters) < 1
                obj.parameters{1} = 0;
            end

            if numel(obj.parameters) < 2
                obj.parameters{2} = 1;
            end
        end
        %% Calculates the Objective Function value
        function fLETxDose = computeLETxDoseObjectiveFunction(obj,LETxDose)
            fLETxDose =  obj.objectiveExpoDiff(LETxDose);
        end
        
        %% Calculates the Objective Function gradient
        function fLETxDoseGrad   = computeLETxDoseObjectiveGradient(obj,LETxDose)
            fLETxDoseGrad = obj.gradientExpoDiff(LETxDose);
        end
    end
     methods (Access = protected)
        function fLETxDose = objectiveExpoDiff(obj,LETxDose)
            % Formula for the Objective Function value
            fLETxDose = (mean(LETxDose(:)) - obj.parameters{1})^obj.parameters{2}; 
        end

        function fLETxDoseGrad = gradientExpoDiff(obj,LETxDose)
            % Derivative of the Objective Function
            fLETxDoseGrad = obj.parameters{2}*((mean(LETxDose(:))-obj.parameters{1}))^(obj.parameters{2}-1) * ones(size(LETxDose(:)))/numel(LETxDose);
        end

    end
    
end


