classdef matRad_SquaredOverdosingLETxDose < LETxDoseObjectives.matRad_LETxDoseObjective
% matRad_SquaredOverdosingLETxDose implements a penalized squared overdosing LETxDose objective
%   See matRad_LETxDoseObjective for interface description
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
        name = 'Squared Overdosing LETxDose';
        parameterNames = {'LETxd^{max}'};
        parameterTypes = {'LETxd'};
    end
    
    properties
        parameters = {30};
        penalty = 1;
    end
    
    methods
        function obj = matRad_SquaredOverdosingLETxDose(penalty,LETxDoseMax)
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
                if nargin == 2 && isscalar(LETxDoseMax)
                    obj.parameters{1} = LETxDoseMax;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
        
        %% Calculates the Objective Function value
        function fLETxDose = computeLETxDoseObjectiveFunction(obj,LETxDose)
            % overLETxDose : LETxDose minus prefered LETxDose
            overLETxDose = LETxDose - obj.parameters{1};
            
            % apply positive operator
            overLETxDose(overLETxDose<0) = 0;
            
            % calculate objective function
            fLETxDose = 1/numel(LETxDose) * (overLETxDose'*overLETxDose);
        end
        
        %% Calculates the Objective Function gradient
        function fLETxDoseGrad   = computeLETxDoseObjectiveGradient(obj,LETxDose)
            % overLETxDose : LETxDose minus prefered LETxDose
            overLETxDose = LETxDose - obj.parameters{1};
            
            % apply positive operator
            overLETxDose(overLETxDose<0) = 0;
            
            % calculate delta
            fLETxDoseGrad = 2 * 1/numel(LETxDose) * overLETxDose;
        end
    end
    
end
