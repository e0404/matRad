classdef matRad_SquaredUnderDADR < DADRObjectives.matRad_DADRObjective
% matRad_SquaredUnderdosing Implements a penalized squared underdosing objective
%   See matRad_DoseObjective for interface description
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
        name = 'Squared UnderDADR';
        parameterNames = {'DADR^{min}'};
        parameterTypes = {'DADR'};
    end
    
    properties
        parameters = {60};
        penalty = 1;
    end
    
    methods
        function obj = matRad_SquaredUnderDADR(penalty,dadrMin)
            %If we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DADRObjectives.matRad_DADRObjective(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                if nargin == 2 && isscalar(dadrMin)
                    obj.parameters{1} = dadrMin;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
        
        %% Calculates the Objective Function value
        function fDADR = computeDADRObjectiveFunction(obj,DADR)
            % overdose : dose minus prefered dose
            underdose = DADR - obj.parameters{1};
            
            % apply positive operator
            underdose(underdose>0) = 0;
            
            % claculate objective function
            fDADR = obj.penalty/numel(DADR) * (underdose'*underdose);
        end
        
        %% Calculates the Objective Function gradient
        function fDADRGrad   = computeDADRObjectiveGradient(obj,DADR)
            % overdose : dose minus prefered dose
            underdose = DADR - obj.parameters{1};
            
            % apply positive operator
            underdose(underdose>0) = 0;
            
            % calculate delta
            fDADRGrad = 2 * obj.penalty/numel(DADR) * underdose;
        end
    end
    
end
