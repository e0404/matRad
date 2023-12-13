classdef matRad_SquaredOverdosingmLETDose < mLETDoseObjectives.matRad_mLETDoseObjective
% matRad_SquaredOverdosingmLETDose Implements a penalized squared overdosing mLETDose objective
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
        name = 'Squared Overdosing mLETDose';
        parameterNames = {'mLETd^{max}'};
        parameterTypes = {'mLETd'};
    end
    
    properties
        parameters = {30};
        penalty = 1;
    end
    
    methods
        function obj = matRad_SquaredOverdosingmLETDose(penalty,mLETDoseMax)
            %If we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@mLETDoseObjectives.matRad_mLETDoseObjective(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                if nargin == 2 && isscalar(mLETDoseMax)
                    obj.parameters{1} = mLETDoseMax;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
        
        %% Calculates the Objective Function value
        function fmLETDose = computemLETDoseObjectiveFunction(obj,mLETDose)
            % overmLETDose : mLETDose minus prefered mLETDose
            overmLETDose = mLETDose - obj.parameters{1};
            
            % apply positive operator
            overmLETDose(overmLETDose<0) = 0;
            
            % calculate objective function
            fmLETDose = 1/numel(mLETDose) * (overmLETDose'*overmLETDose);
        end
        
        %% Calculates the Objective Function gradient
        function fmLETDoseGrad   = computemLETDoseObjectiveGradient(obj,mLETDose)
            % overmLETDose : mLETDose minus prefered mLETDose
            overmLETDose = mLETDose - obj.parameters{1};
            
            % apply positive operator
            overmLETDose(overmLETDose<0) = 0;
            
            % calculate delta
            fmLETDoseGrad = 2 * 1/numel(mLETDose) * overmLETDose;
        end
    end
    
end
