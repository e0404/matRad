classdef matRad_MeanDirtyDose < DirtyDoseObjectives.matRad_DirtyDoseObjective
% matRad_MeanDirtyDose Implements a penalized MeanDirtyDose objective
%   See matRad_DirtyDoseObjective for interface description
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
        name = 'Mean Dirty Dose';
        parameterNames = {'d^{ref}','f_{diff}'}; %When optimizing to a reference, one might consider using a quadratic relationship with a non-linear optimizer
        parameterTypes = {'dirtyDose',{'Linear','Quadratic'}};
    end
    
    properties
        parameters = {0,1};        
        penalty = 1;
    end
    
    methods 
        function obj = matRad_MeanDirtyDose(penalty,dMeanRef,fDiff)
           
            % if we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DirtyDoseObjectives.matRad_DirtyDoseObjective(inputStruct);
            
            if ~initFromStruct
                if nargin < 3 || ~ischar(fDiff)
                    fDiff = 'Linear';
                end
                
                fDiffIx = find(strcmp(fDiff,obj.parameterTypes{2}));
                
                if isempty(fDiffIx) || numel(fDiffIx) > 1
                    fDiffIx = 1;
                    matRad_cfg = MatRad_Config.instance();                    
                    matRad_cfg.dispWarning('Mean dirty dose difference function can only be %s! Using %s difference.', strjoin(obj.parameterTypes{2},' or '), obj.parameterTypes{2}{fDiffIx});
                end
                
                obj.parameters{2} = fDiffIx;


                if nargin >= 2 && isscalar(dMeanRef)
                    obj.parameters{1} = dMeanRef;
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
        function fDirtyDose = computeDirtyDoseObjectiveFunction(obj,dirtyDose)
            switch obj.parameters{2}
                case 1
                    fDirtyDose =  obj.objectiveLinearDiff(dirtyDose);
                case 2
                    fDirtyDose =  obj.objectiveQuadraticDiff(dirtyDose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid setting for %s in Mean Dirty Dose Objective!',obj.parameterNames{2});  
            end
        end
        
        %% Calculates the Objective Function gradient
        function fDirtyDoseGrad   = computeDirtyDoseObjectiveGradient(obj,dirtyDose)
            switch obj.parameters{2}
                case 1
                    fDirtyDoseGrad = obj.gradientLinearDiff(dirtyDose);
                case 2
                    fDirtyDoseGrad = obj.gradientQuadraticDiff(dirtyDose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid setting for %s in Mean Dirty Dose Objective!',obj.parameterNames{2});  
            end
        end
    end

    methods (Access = protected)
        function fDirtyDose = objectiveQuadraticDiff(obj,dirtyDose)
            fDirtyDose = (mean(dirtyDose(:)) - obj.parameters{1})^2;
        end

        function fDirtyDoseGrad = gradientQuadraticDiff(obj,dirtyDose)
            fDirtyDoseGrad = 2*(mean(dirtyDose(:))-obj.parameters{1}) * ones(size(dirtyDose(:)))/numel(dirtyDose);
        end

        function fDirtyDose = objectiveLinearDiff(obj,dirtyDose)
            fDirtyDose = abs(mean(dirtyDose(:)) - obj.parameters{1});
        end

        function fDirtyDoseGrad = gradientLinearDiff(obj,dirtyDose)
            fDirtyDoseGrad = (1/numel(dirtyDose))*sign(dirtyDose(:)-obj.parameters{1});
        end
    end
    
end

