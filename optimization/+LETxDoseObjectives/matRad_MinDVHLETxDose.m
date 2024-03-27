classdef matRad_MinDVHLETxDose < LETxDoseObjectives.matRad_LETxDoseObjective
% matRad_MinDVHLETxDose Implements a penalized MinDVHLETxDose objective
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
        name = 'Min DVH LETxDose';
        parameterNames = {'LETxDose', 'V^{min}'};
        parameterTypes = {'LETxDose','numeric'};
    end
    
    properties
        parameters = {60,95};
        penalty = 1;
    end
    
    methods
        function obj = matRad_MinDVHLETxDose(penalty,LETxDoseRef,vMinPercent)
            
            % if we have a struct in first argument
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
                if nargin >= 3 && isscalar(vMinPercent)
                    obj.parameters{2} = vMinPercent;
                end
                
                if nargin >= 2 && isscalar(LETxDoseRef)
                    obj.parameters{1} = LETxDoseRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
            
        end        
        %% Calculates the Objective Function value
        function fLETxDose = computeLETxDoseObjectiveFunction(obj,LETxDose)                       
            % get reference Volume
            refVol = obj.parameters{2}/100;
            
            % calc deviation
            deviation = LETxDose - obj.parameters{1};

            % calc d_ref2: V(d_ref2) = refVol
            d_ref2 = matRad_calcInversDVHLETxDose(refVol,LETxDose);

            
            deviation(LETxDose > obj.parameters{1} | LETxDose < d_ref2) = 0;
   
            % calculate objective function
            fLETxDose = (1/numel(LETxDose))*(deviation'*deviation);
        end
        
        %% Calculates the Objective Function gradient
        function fLETxDoseGrad   = computeLETxDoseObjectiveGradient(obj,LETxDose)
            % get reference Volume
            refVol = obj.parameters{2}/100;
            
            % calc deviation
            deviation = LETxDose - obj.parameters{1};
            
            % calc d_ref2: V(d_ref2) = refVol
            d_ref2 = matRad_calcInversDVHLETxDose(refVol,LETxDose);
            
            deviation(LETxDose > obj.parameters{1} | LETxDose < d_ref2) = 0;

            % calculate delta
            fLETxDoseGrad = (2/numel(LETxDose))*deviation;
        end
    end
    
end