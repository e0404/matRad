classdef matRad_MaxDVHLETd < LETdObjectives.matRad_LETdObjective
% matRad_MaxDVHLETd Implements a penalized maximum DVH LETd objective
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
        name = 'Max DVH LETd';
        parameterNames = {'LETd', 'V^{max}'};
        parameterTypes = {'LETd','numeric'};
    end
    
    properties
        parameters = {30,95};
        penalty = 1;
    end
    
    methods 
        function obj = matRad_MaxDVHLETd(penalty,LETdRef,vMaxPercent)
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
                if nargin >= 3 && isscalar(vMaxPercent)
                    obj.parameters{2} = vMaxPercent;
                end
                
                if nargin >= 2 && isscalar(LETdRef)
                    obj.parameters{1} = LETdRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end        
        
        %% Calculates the Objective Function value
        function fLETd = computeLETdObjectiveFunction(obj,LETd)                       
            % get reference Volume
            refVol = obj.parameters{2}/100;
            
            % calc deviation
            deviation = LETd - obj.parameters{1};

            % calc d_ref2: V(d_ref2) = refVol
            d_ref2 = matRad_calcInversDVHLET(refVol,LETd);

            
            deviation(LETd < obj.parameters{1} | LETd > d_ref2) = 0;
   
            % calculate objective function
            fLETd = (1/numel(LETd))*(deviation'*deviation);
        end
        
        %% Calculates the Objective Function gradient
        function fLETdGrad   = computeLETdObjectiveGradient(obj,LETd)
            % get reference Volume
            refVol = obj.parameters{2}/100;
            
            % calc deviation
            deviation = LETd - obj.parameters{1};
            
            % calc d_ref2: V(d_ref2) = refVol
            d_ref2 = matRad_calcInversDVHLETd(refVol,LETd);
            
            deviation(LETd < obj.parameters{1} | LETd > d_ref2) = 0;

            % calculate delta
            fLETdGrad = (2/numel(LETd))*deviation;
        end
    end
    
end
