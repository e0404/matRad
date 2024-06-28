classdef matRad_SquaredDeviationDirtyDose < DirtyDoseObjectives.matRad_DirtyDoseObjective
% matRad_SquaredDeviationDirtyDose Implements a penalized least squares dirtyDose objective
%   See matRad_DirtyDoseObjective for interface description
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
        name = 'Squared Deviation Dirty Dose';
        parameterNames = {'d^{ref}'};
        parameterTypes = {'dirtyDose'};
    end
    
    properties
        parameters = {60};
        penalty = 1;
    end

    
    methods
        function obj = matRad_SquaredDeviationDirtyDose(penalty,dRef)
            
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
        function fDirtyDose = computeDirtyDoseObjectiveFunction(obj,dirtyDose)
            % deviation : dirtyDose minus prefered dose
            deviation = dirtyDose - obj.parameters{1};

            % calculate objective function
            fDirtyDose = 1/numel(dirtyDose) * (deviation'*deviation);
        end
        
        %% Calculates the Objective Function gradient
        function fDirtyDoseGrad   = computeDirtyDoseObjectiveGradient(obj,dirtyDose)
            % deviation : dirtyDose minus prefered dose
            deviation = dirtyDose - obj.parameters{1};
            
            % calculate delta
            fDirtyDoseGrad = 2 * 1/numel(dirtyDose) * deviation;
        end
    end
end