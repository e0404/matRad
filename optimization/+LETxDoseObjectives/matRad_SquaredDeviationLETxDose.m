classdef matRad_SquaredDeviationLETxDose < LETxDoseObjectives.matRad_LETxDoseObjective
% matRad_SquaredDeviationLETxDose implements a penalized least squares objective
%   See matRad_LETxDoseObjective for interface description
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
        name = 'Squared Deviation LETxDose';
        parameterNames = {'LETxd^{ref}'};
        parameterTypes = {'LETxd'};
    end
    
    properties
        parameters = {60};
        penalty = 1;
    end
    
    methods
        function obj = matRad_SquaredDeviationLETxDose(penalty,LETxDoseRef)
            
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
            % deviation : mLETDose minus prefered mLETDose
            deviation = LETxDose - obj.parameters{1};
            % calculate objective function
            fLETxDose = 1/numel(LETxDose) * (deviation'*deviation);
        end
        
        %% Calculates the Objective Function gradient
        function fLETxDoseGrad   = computeLETxDoseObjectiveGradient(obj,LETxDose)
            % deviation : LETxDose minus prefered LETxDose
            deviation = LETxDose - obj.parameters{1};
            
            % calculate delta
            fLETxDoseGrad = 2 * 1/numel(LETxDose) * deviation;
        end
    end
    
end