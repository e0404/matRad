classdef matRad_ExpoMeanDose < DoseObjectives.matRad_DoseObjective
% matRad_ExpoMeanDose implements a penalized exponential mean dose objective
% f = ((sum(d)/N)-dRef)^k
%   See matRad_DoseObjective for interface description
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
        name = 'Exponential Mean Dose';
        parameterNames = {'d^{ref}','k'};       % k is the exponent of the mean function
        parameterTypes = {'dose','numeric'};
    end
    
    properties
        parameters = {0,2};        
        penalty = 1;
    end
    
    methods 
        function obj = matRad_ExpoMeanDose(penalty,dMeanRef,k)
           
            % if we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            % Call Superclass Constructor (for struct initialization)
            obj@DoseObjectives.matRad_DoseObjective(inputStruct);
            
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
        function fDose = computeDoseObjectiveFunction(obj,dose)
            fDose =  obj.objectiveExpoDiff(dose);
        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad   = computeDoseObjectiveGradient(obj,dose)
            fDoseGrad = obj.gradientExpoDiff(dose);
        end
    end
     methods (Access = protected)
        function fDose = objectiveExpoDiff(obj,dose)
            % Formula for the Objective Function value
            fDose = (mean(dose(:)) - obj.parameters{1})^obj.parameters{2}; 
        end

        function fDoseGrad = gradientExpoDiff(obj,dose)
            % Derivative of the Objective Function
            fDoseGrad = obj.parameters{2}*((mean(dose(:))-obj.parameters{1}))^(obj.parameters{2}-1) * ones(size(dose(:)))/numel(dose);
        end

    end
    
end


