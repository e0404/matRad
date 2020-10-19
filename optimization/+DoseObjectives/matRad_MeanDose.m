classdef matRad_MeanDose < DoseObjectives.matRad_DoseObjective
% matRad_MeanDose Implements a penalized MeanDose objective
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
        name = 'Mean Dose';
        parameterNames = {'d^{ref}'};
        parameterTypes = {'dose'};
        %parameterNames = {};
        %parameterIsDose = [];
    end
    
    properties
        parameters = {0};        
        penalty = 1;
    end
    
    methods 
        function obj = matRad_MeanDose(penalty,dMeanRef)
           
            % if we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DoseObjectives.matRad_DoseObjective(inputStruct);
            
            if ~initFromStruct
                if nargin == 2 && isscalar(dMeanRef)
                    obj.parameters{1} = dMeanRef;
                end

                if nargin >= 1 && isscalar(s)
                    obj.penalty = penalty;
                end
            end
            
        end       
        
        %% Calculates the Objective Function value
        function fDose = computeDoseObjectiveFunction(obj,dose)
            %fDose = obj.penalty * abs(mean(dose(:)) - obj.parameters{1});
            fDose = obj.penalty * (mean(dose(:)) - obj.parameters{1})^2;
        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad   = computeDoseObjectiveGradient(obj,dose)
            %fDoseGrad = (obj.penalty/numel(dose))*sign(dose(:)-obj.parameters{1});
            fDoseGrad = obj.penalty*2*(mean(dose(:))-obj.parameters{1}) * ones(size(dose(:)))/numel(dose);
        end
    end
    
end

