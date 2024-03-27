classdef matRad_SquaredOverdosingLETd < LETdObjectives.matRad_LETdObjective
% matRad_SquaredOverdosingLETd Implements a penalized squared overdosing LETd objective
%   See matRad_LETdObjective for interface description
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
        name = 'Squared Overdosing LETd';
        parameterNames = {'LETd^{max}'};
        parameterTypes = {'LETd'};
    end
    
    properties
        parameters = {30};
        penalty = 1;
    end
    
    methods
        function obj = matRad_SquaredOverdosingLETd(penalty,LETdMax)
            %If we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@LETdObjectives.matRad_LETdObjective(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                if nargin == 2 && isscalar(LETdMax)
                    obj.parameters{1} = LETdMax;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
        
        %% Calculates the Objective Function value
        function fLETd = computeLETdObjectiveFunction(obj,LETd)
            % overLETd : LETd minus prefered dose
            overLETd = LETd - obj.parameters{1};
            
            % apply positive operator
            overLETd(overLETd<0) = 0;
            
            % calculate objective function
            fLETd = 1/numel(LETd) * (overLETd'*overLETd);
        end
        
        %% Calculates the Objective Function gradient
        function fLETdGrad   = computeLETdObjectiveGradient(obj,LETd)
            % overLETd : LETd minus prefered dose
            overLETd = LETd - obj.parameters{1};
            
            % apply positive operator
            overLETd(overLETd<0) = 0;
            
            % calculate delta
            fLETdGrad = 2 * 1/numel(LETd) * overLETd;
        end
    end
    
end
