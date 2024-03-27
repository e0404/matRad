classdef matRad_LETxDoseVariance < LETxDoseObjectives.matRad_LETxDoseObjective
% matRad_LETxDoseVariance Implements a variance objective for homogenous LETxDose
% See matRad_LETxDoseObjective for interface description
%
% References 
%     -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
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
        name = 'Variance';
        parameterNames = {};
        parameterTypes = {};
    end
    
    properties
        parameters = {};
        penalty = 1;
    end
    
    methods
        function obj = matRad_LETxDoseVariance(penalty)
            
            %If we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@LETxDoseObjectives.matRad_LETxDoseObjective(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                if nargin == 2 && isscalar(LETxDoseRef)
                    obj.parameters{1} = LETxDoseRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
   
        %% Calculates the Objective Function value
        function fLETxDose = computeLETxDoseObjectiveFunction(obj,LETxDose)
                 fLETxDose = obj.penalty * var(LETxDose);
            
        end
        
        %% Calculates the Objective Function gradient
        function fLETxDoseGrad   = computeLETxDoseObjectiveGradient(obj,LETxDose)
            
            fLETxDoseGrad = obj.penalty * 2/(numel(LETxDose) - 1) * (LETxDose - mean(LETxDose));
        end
    end
    
end