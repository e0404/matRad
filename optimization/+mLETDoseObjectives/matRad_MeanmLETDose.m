classdef matRad_MeanmLETDose < mLETDoseObjectives.matRad_mLETDoseObjective
% matRad_MeanmLETDose Implements a penalized MeanmLETDose objective
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
        name = 'Mean mLETDose';
        parameterNames = {'mLETDose^{ref}','f_{diff}'}; %When optimizing to a reference, one might consider using a quadratic relationship with a non-linear optimizer
        parameterTypes = {'mLETDose',{'Linear','Quadratic'}};
    end
    
    properties
        parameters = {0,1};        
        penalty = 1;
    end
    
    methods 
        function obj = matRad_MeanmLETDose(penalty,mLETDoseRef,fDiff)
           
            % if we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@mLETDoseObjectives.matRad_mLETDoseObjective(inputStruct);
            
            if ~initFromStruct
                if nargin < 3 || ~ischar(fDiff)
                    fDiff = 'Linear';
                end
                
                fDiffIx = find(strcmp(fDiff,obj.parameterTypes{2}));
                
                if isempty(fDiffIx) || numel(fDiffIx) > 1
                    fDiffIx = 1;
                    matRad_cfg = MatRad_Config.instance();                    
                    matRad_cfg.dispWarning('Mean mLETDose difference function can only be %s! Using %s difference.', strjoin(obj.parameterTypes{2},' or '), obj.parameterTypes{2}{fDiffIx});
                end
                
                obj.parameters{2} = fDiffIx;


                if nargin >= 2 && isscalar(mLETDoseRef)
                    obj.parameters{1} = mLETDoseRef;
                end

                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        
        %% Downwards compatability / set default values
            %TODO: maybe move into set method for parameters
            if numel(obj.parameters) < 1
                obj.parameters{1} = 0;
            end

            if numel(obj.parameters) < 2
                obj.parameters{2} = 1;
            end
            
        end       
        
        %% Calculates the Objective Function value
        function fmLETDose = computemLETDoseObjectiveFunction(obj,mLETDose)
            switch obj.parameters{2}
                case 1
                    fmLETDose =  obj.objectiveLinearDiff(mLETDose);
                case 2
                    fmLETDose =  obj.objectiveQuadraticDiff(mLETDose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid setting for %s in Mean mLETDose Objective!',obj.parameterNames{2});  
            end
        end
        
        %% Calculates the Objective Function gradient
        function fmLETDoseGrad   = computemLETDoseObjectiveGradient(obj,mLETDose)
            switch obj.parameters{2}
                case 1
                    fmLETDoseGrad = obj.gradientLinearDiff(mLETDose);
                case 2
                    fmLETDoseGrad = obj.gradientQuadraticDiff(mLETDose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid setting for %s in Mean mLETDose Objective!',obj.parameterNames{2});  
            end
        end
    end

    methods (Access = protected)
        function fmLETDose = objectiveQuadraticDiff(obj,mLETDose)
            fmLETDose = (mean(mLETDose(:)) - obj.parameters{1})^2;
        end

        function fmLETDoseGrad = gradientQuadraticDiff(obj,mLETDose)
            fmLETDoseGrad = 2*(mean(mLETDose(:))-obj.parameters{1}) * ones(size(mLETDose(:)))/numel(mLETDose);
        end

        function fmLETDose = objectiveLinearDiff(obj,mLETDose)
            fmLETDose = abs(mean(mLETDose(:)) - obj.parameters{1});
        end

        function fmLETDoseGrad = gradientLinearDiff(obj,mLETDose)
            fmLETDoseGrad = (1/numel(mLETDose))*sign(mLETDose(:)-obj.parameters{1});
        end
    end
    
end

