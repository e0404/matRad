classdef matRad_SquaredDeviationDADR < DADRObjectives.matRad_DADRObjective
% matRad_SquaredDeviation Implements a penalized least squares objective
%   See matRad_DoseObjective for interface description
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
        name = 'Squared Deviation';
        parameterNames = {'DADR^{ref}'};
        parameterTypes = {'DADR'};
    end
    
    properties
        parameters = {60};
        penalty = 1;
    end
    
    methods
        function obj = matRad_SquaredDeviationDADR(penalty,dadrRef)
            
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
                if nargin == 2 && isscalar(dadrRef)
                    obj.parameters{1} = dadrRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
   
        %% Calculates the Objective Function value
        function fDADR = computeDADRObjectiveFunction(obj,DADR)
            % deviation : dose minus prefered dose
            deviation = DADR - obj.parameters{1};
            % claculate objective function
            fDADR = obj.penalty/numel(DADR) * (deviation'*deviation);
        end
        
        %% Calculates the Objective Function gradient
        function fDADRGrad   = computeDADRObjectiveGradient(obj,DADR)
            % deviation : Dose minus prefered dose
            deviation = DADR - obj.parameters{1};
            
            % calculate delta
            fDADRGrad = 2 * obj.penalty/numel(DADR) * deviation;
        end
    end
    
end
