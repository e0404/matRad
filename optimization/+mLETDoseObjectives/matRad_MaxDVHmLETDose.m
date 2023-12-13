classdef matRad_MaxDVHmLETDose < mLETDoseObjectives.matRad_mLETDoseObjective
% matRad_MaxDVHmLETDose Implements a penalized maximum DVH mLETDose objective
%   See matRad_mLETDoseObjective for interface description
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
        name = 'Max DVH mLETDose';
        parameterNames = {'mLETDose', 'V^{max}'};
        parameterTypes = {'mLETDose','numeric'};
    end
    
    properties
        parameters = {30,95};
        penalty = 1;
    end
    
    methods 
        function obj = matRad_MaxDVHmLETDose(penalty,mLETDoseRef,vMaxPercent)
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
                if nargin >= 3 && isscalar(vMaxPercent)
                    obj.parameters{2} = vMaxPercent;
                end
                
                if nargin >= 2 && isscalar(mLETDoseRef)
                    obj.parameters{1} = mLETDoseRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end        
        
        %% Calculates the Objective Function value
        function fmLETDose = computemLETDoseObjectiveFunction(obj,mLETDose)                       
            % get reference Volume
            refVol = obj.parameters{2}/100;
            
            % calc deviation
            deviation = mLETDose - obj.parameters{1};

            % calc d_ref2: V(d_ref2) = refVol
            d_ref2 = matRad_calcInversDVHmLETDose(refVol,mLETDose);

            
            deviation(mLETDose < obj.parameters{1} | mLETDose > d_ref2) = 0;
   
            % claculate objective function
            fmLETDose = (1/numel(mLETDose))*(deviation'*deviation);
        end
        
        %% Calculates the Objective Function gradient
        function fmLETDoseGrad   = computemLETDoseObjectiveGradient(obj,mLETDose)
            % get reference Volume
            refVol = obj.parameters{2}/100;
            
            % calc deviation
            deviation = mLETDose - obj.parameters{1};
            
            % calc d_ref2: V(d_ref2) = refVol
            d_ref2 = matRad_calcInversDVHmLETDose(refVol,mLETDose);
            
            deviation(mLETDose < obj.parameters{1} | mLETDose > d_ref2) = 0;

            % calculate delta
            fmLETDoseGrad = (2/numel(mLETDose))*deviation;
        end
    end
    
end
