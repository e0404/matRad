classdef matRad_DirtyDoseVariance < DirtyDoseObjectives.matRad_DirtyDoseObjective
% matRad_DirtyDoseVariance implements a variance objective for
% homogenous dirty dose distribution
%   See matRad_DoseObjective for interface description
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
        function obj = matRad_DirtyDoseVariance(penalty)
            
            %If we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DirtyDoseObjectives.matRad_DirtyDoseObjective(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                if nargin == 2 && isscalar(dRef)
                    obj.parameters{1} = dRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
   
        %% Calculates the Objective Function value
        function fDirtyDose       = computeDirtyDoseObjectiveFunction(obj,dirtyDose)
                 fDirtyDose = obj.penalty * var(dirtyDose);
            
        end
        
        %% Calculates the Objective Function gradient
        function fDirtyDoseGrad   = computeDirtyDoseObjectiveGradient(obj,dirtyDose)
                 fDirtyDoseGrad = obj.penalty * 2/(numel(dirtyDose) - 1) * (dirtyDose - mean(dirtyDose));
        end
    end
    
end