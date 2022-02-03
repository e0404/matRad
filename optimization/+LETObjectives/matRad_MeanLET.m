classdef matRad_MeanLET < LETObjectives.matRad_LETObjective
% matRad_MeanLET Implements a penalized MeanLET objective
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
        name = 'Mean LET';
        %parameterNames = {'LET^{ref}'};
        %parameterTypes = {'LET'};
        parameterNames = {};
        parameterTypes = {};
        %parameterIsLET = [];
    end
    
    properties
        parameters = {};        
        penalty = 1;
    end
    
    methods 
        %function obj = matRad_MeanLET(penalty,LETMeanRef)
        function obj = matRad_MeanLET(penalty,LETMeanRef)
           
            % if we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@LETObjectives.matRad_LETObjective(inputStruct);
            
            if ~initFromStruct
                if nargin >= 1 && isscalar(inputStruct)
                    obj.penalty = penalty;
                end
            end
            
        end       
        
        %% Calculates the Objective Function value
        function fLET = computeLETObjectiveFunction(obj,LET)
            %fLET = obj.penalty * abs(mean(LET(:)) - obj.parameters{1});
            fLET = obj.penalty * mean(LET(:));
        end
        
        %% Calculates the Objective Function gradient
        function fLETGrad   = computeLETObjectiveGradient(obj,LET)
            %fLETGrad = (obj.penalty/numel(LET))*sign(LET(:)-obj.parameters{1});
            fLETGrad = obj.penalty * ones(size(LET))./numel(LET);
        end
    end
    
end

