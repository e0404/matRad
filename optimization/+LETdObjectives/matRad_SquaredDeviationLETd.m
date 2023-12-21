classdef matRad_SquaredDeviationLETd < LETdObjectives.matRad_LETdObjective
% matRad_SquaredDeviationLETd Implements a penalized least squares objective
%   See matRad_LETdObjective for interface description
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
        name = 'Squared Deviation LETd';
        parameterNames = {'LETd^{ref}'};
        parameterTypes = {'LETd'};
    end
    
    properties
        parameters = {60};
        penalty = 1;
    end

    
    methods
        function obj = matRad_SquaredDeviationLETd(penalty,LETdRef)
            
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
                if nargin == 2 && isscalar(LETdRef)
                    obj.parameters{1} = LETdRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
   
        %% Calculates the Objective Function value
        function fLETd = computeLETdObjectiveFunction(obj,LETd)
            % deviation : LETd minus prefered LETd
            deviation = LETd - obj.parameters{1};
            % calculate objective function
            fLETd = 1/numel(LETd) * (deviation'*deviation);
        end
        
        %% Calculates the Objective Function gradient
        function fLETdGrad   = computeLETdObjectiveGradient(obj,LETd)
            % deviation : LETd minus prefered LETd
            deviation = LETd - obj.parameters{1};
            
            % calculate delta
            fLETdGrad = 2/numel(LETd) * deviation;
        end
    end
    
end
