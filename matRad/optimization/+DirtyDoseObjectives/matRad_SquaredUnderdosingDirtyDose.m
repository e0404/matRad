classdef matRad_SquaredUnderdosingDirtyDose < DirtyDoseObjectives.matRad_DirtyDoseObjective
% matRad_SquaredUnderdosingDirtyDose implements a penalized dirty dose
%   See matRad_DirtyDoseObjective for interface description
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
    
    properties (Constant)
        name = 'Squared Underdosing Dirty Dose';
        parameterNames = {'d^{min}'};
        parameterTypes = {'dirtyDose'};
    end
    
    properties
        parameters = {60};
        penalty = 1;
    end
    
    methods
        function obj = matRad_SquaredUnderdosingDirtyDose(penalty,dMin)
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
                if nargin >= 2 && isscalar(dMin)
                    obj.parameters{1} = dMin;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
        
        %% Calculates the Objective Function value
        function fDirtyDose = computeDirtyDoseObjectiveFunction(obj,dirtyDose)
            % underdose : dirtyDose minus prefered dose
            underdose = dirtyDose - obj.parameters{1};
            
            % apply positive operator
            underdose(underdose>0) = 0;
            
            % calculate objective function
            fDirtyDose = 1/numel(dirtyDose) * (underdose'*underdose);
        end
        
        %% Calculates the Objective Function gradient
        function fDirtyDoseGrad   = computeDirtyDoseObjectiveGradient(obj,dirtyDose)
            % underdose : dirtyDose minus prefered dose
            underdose = dirtyDose - obj.parameters{1};
            
            % apply positive operator
            underdose(underdose>0) = 0;
            
            % calculate delta
            fDirtyDoseGrad = 2/numel(dirtyDose) * underdose;
        end
    end
    
end
