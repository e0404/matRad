classdef matRad_SquaredUnderLET < LETObjectives.matRad_LETObjective
% matRad_SquaredUnderLETing Implements a penalized squared underLETing objective
%   See matRad_LETObjective for interface description
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
        name = 'Squared UnderLET';
        parameterNames = {'LETd^{min}'};
        parameterTypes = {'LETd'};
    end
    
    properties
        parameters = {60};
        penalty = 1;
    end
    
    methods
        function obj = matRad_SquaredUnderLET(penalty,LETMin)
            %If we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@LETObjectives.matRad_LETObjective(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                if nargin == 2 && isscalar(LETMin)
                    obj.parameters{1} = LETMin;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
        
        %% Calculates the Objective Function value
        function fLET = computeLETObjectiveFunction(obj,LET)
            % overLET : LET minus prefered LET
            underLET = LET - obj.parameters{1};
            
            % apply positive operator
            underLET(underLET>0) = 0;
            
            % claculate objective function
            fLET = obj.penalty/numel(LET) * (underLET'*underLET);
        end
        
        %% Calculates the Objective Function gradient
        function fLETGrad   = computeLETObjectiveGradient(obj,LET)
            % overLET : LET minus prefered LET
            underLET = LET - obj.parameters{1};
            
            % apply positive operator
            underLET(underLET>0) = 0;
            
            % calculate delta
            fLETGrad = 2 * obj.penalty/numel(LET) * underLET;
        end
    end
    
end
