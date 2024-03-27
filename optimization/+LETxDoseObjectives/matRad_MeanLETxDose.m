classdef matRad_MeanLETxDose < LETxDoseObjectives.matRad_LETxDoseObjective
% matRad_MeanLETxDose implements a penalized MeanLETxDose objective
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
        name = 'Mean LETxDose';
        parameterNames = {'LETxDose^{ref}','f_{diff}'}; %When optimizing to a reference, one might consider using a quadratic relationship with a non-linear optimizer
        parameterTypes = {'LETxDose',{'Linear','Quadratic'}};
    end
    
    properties
        parameters = {0,1};        
        penalty = 1;
    end
    
    methods 
        function obj = matRad_MeanLETxDose(penalty,LETxDoseRef,fDiff)
           
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
            
            if ~initFromStruct
                if nargin < 3 || ~ischar(fDiff)
                    fDiff = 'Linear';
                end
                
                fDiffIx = find(strcmp(fDiff,obj.parameterTypes{2}));
                
                if isempty(fDiffIx) || numel(fDiffIx) > 1
                    fDiffIx = 1;
                    matRad_cfg = MatRad_Config.instance();                    
                    matRad_cfg.dispWarning('Mean LETxDose difference function can only be %s! Using %s difference.', strjoin(obj.parameterTypes{2},' or '), obj.parameterTypes{2}{fDiffIx});
                end
                
                obj.parameters{2} = fDiffIx;


                if nargin >= 2 && isscalar(LETxDoseRef)
                    obj.parameters{1} = LETxDoseRef;
                end

                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        
        %% Downwards compatability / set default values
            
            if numel(obj.parameters) < 1
                obj.parameters{1} = 0;
            end

            if numel(obj.parameters) < 2
                obj.parameters{2} = 1;
            end
            
        end       
        
        %% Calculates the Objective Function value
        function fLETxDose = computeLETxDoseObjectiveFunction(obj,LETxDose)
            switch obj.parameters{2}
                case 1
                    fLETxDose =  obj.objectiveLinearDiff(LETxDose);
                case 2
                    fLETxDose =  obj.objectiveQuadraticDiff(LETxDose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid setting for %s in Mean LETxDose Objective!',obj.parameterNames{2});  
            end
        end
        
        %% Calculates the Objective Function gradient
        function fLETxDoseGrad   = computeLETxDoseObjectiveGradient(obj,LETxDose)
            switch obj.parameters{2}
                case 1
                    fLETxDoseGrad = obj.gradientLinearDiff(LETxDose);
                case 2
                    fLETxDoseGrad = obj.gradientQuadraticDiff(LETxDose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid setting for %s in Mean LETxDose Objective!',obj.parameterNames{2});  
            end
        end
    end

    methods (Access = protected)
        function fLETxDose = objectiveQuadraticDiff(obj,LETxDose)
            fLETxDose = (mean(LETxDose(:)) - obj.parameters{1})^2;
        end

        function fLETxDoseGrad = gradientQuadraticDiff(obj,LETxDose)
            fLETxDoseGrad = 2*(mean(LETxDose(:))-obj.parameters{1}) * ones(size(LETxDose(:)))/numel(LETxDose);
        end

        function fLETxDose = objectiveLinearDiff(obj,LETxDose)
            fLETxDose = abs(mean(LETxDose(:)) - obj.parameters{1});
        end

        function fLETxDoseGrad = gradientLinearDiff(obj,LETxDose)
            fLETxDoseGrad = (1/numel(mLETxDose))*sign(LETxDose(:)-obj.parameters{1});
        end
    end
    
end

