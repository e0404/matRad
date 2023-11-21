classdef matRad_ExpoMeanDirtyDose < DirtyDoseObjectives.matRad_DirtyDoseObjective
% matRad_ExpoMeanDirtyDose implements a penalized exponential mean dirty dose objective
% f = ((sum(d)/N)-dRef)^k
%   See matRad_DirtyDoseObjective for interface description
%
% References
%   matRad_MeanDirtyDose
%   matRad_SquaredDeviationDirtyDose
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
        name = 'Exponential Mean Dirty Dose';
        parameterNames = {'d^{ref}','k'};       % k is the exponent of the mean function
        parameterTypes = {'dirtyDose','numeric'};
    end
    
    properties
        parameters = {0,2};        
        penalty = 1;
    end
    
    methods 
        function obj = matRad_ExpoMeanDirtyDose(penalty,dMeanRef,k)
           
            % if we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            % Call Superclass Constructor (for struct initialization)
            obj@DirtyDoseObjectives.matRad_DirtyDoseObjective(inputStruct);
            
            % now handle initialization from other parameters
            if ~initFromStruct
                if nargin == 3 && isscalar(k)
                    obj.parameters{2} = k;
                end
                if nargin >= 2 && isscalar(dMeanRef)
                    obj.parameters{1} = dMeanRef;
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
        function fDirtyDose = computeDirtyDoseObjectiveFunction(obj,dirtyDose)
            fDirtyDose =  obj.objectiveExpoDiff(dirtyDose);
        end
        
        %% Calculates the Objective Function gradient
        function fDirtyDoseGrad   = computeDirtyDoseObjectiveGradient(obj,dirtyDose)
            fDirtyDoseGrad = obj.gradientExpoDiff(dirtyDose);
        end
    end
     methods (Access = protected)
        function fDirtyDose = objectiveExpoDiff(obj,dirtyDose)
            % Formula for the Objective Function value
            fDirtyDose = (mean(dirtyDose(:)) - obj.parameters{1})^obj.parameters{2}; 
        end

        function fDirtyDoseGrad = gradientExpoDiff(obj,dirtyDose)
            % Derivative of the Objective Function
            fDirtyDoseGrad = obj.parameters{2}*((mean(dirtyDose(:))-obj.parameters{1}))^(obj.parameters{2}-1) * ones(size(dirtyDose(:)))/numel(dirtyDose);
        end

    end
    
end


