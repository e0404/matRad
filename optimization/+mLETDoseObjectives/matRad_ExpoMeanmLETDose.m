classdef matRad_ExpoMeanmLETDose < mLETDoseObjectives.matRad_mLETDoseObjective
% matRad_ExpoMeanmLETDose implements a penalized exponential mean mLETDose objective
% f = ((sum(mLETDose)/N)-mLETDoseRef)^k
%   See matRad_mLETDoseObjective for interface description
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
        name = 'Exponential Mean mLETDose';
        parameterNames = {'mLETDose^{ref}','k'};       % k is the exponent of the mean function
        parameterTypes = {'mLETDose','numeric'};
    end
    
    properties
        parameters = {0,2};        
        penalty = 1;
    end
    
    methods 
        function obj = matRad_ExpoMeanmLETDose(penalty,mLETDoseRef,k)
           
            % if we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            % Call Superclass Constructor (for struct initialization)
            obj@mLETDoseObjectives.matRad_mLETDoseObjective(inputStruct);
            
            % now handle initialization from other parameters
            if ~initFromStruct
                if nargin == 3 && isscalar(k)
                    obj.parameters{2} = k;
                end
                if nargin >= 2 && isscalar(mLETDoseRef)
                    obj.parameters{1} = mLETDoseRef;
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
        function fmLETDose = computemLETDoseObjectiveFunction(obj,mLETDose)
            fmLETDose =  obj.objectiveExpoDiff(mLETDose);
        end
        
        %% Calculates the Objective Function gradient
        function fmLETDoseGrad   = computemLETDoseObjectiveGradient(obj,mLETDose)
            fmLETDoseGrad = obj.gradientExpoDiff(mLETDose);
        end
    end
     methods (Access = protected)
        function fmLETDose = objectiveExpoDiff(obj,mLETDose)
            % Formula for the Objective Function value
            fmLETDose = (mean(mLETDose(:)) - obj.parameters{1})^obj.parameters{2}; 
        end

        function fmLETDoseGrad = gradientExpoDiff(obj,mLETDose)
            % Derivative of the Objective Function
            fmLETDoseGrad = obj.parameters{2}*((mean(mLETDose(:))-obj.parameters{1}))^(obj.parameters{2}-1) * ones(size(mLETDose(:)))/numel(mLETDose);
        end

    end
    
end


