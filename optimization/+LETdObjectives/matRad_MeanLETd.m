classdef matRad_MeanLETd < LETdObjectives.matRad_LETdObjective
% matRad_MeanLETd Implements a penalized MeanLETd objective
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
        name = 'Mean LETd';
        parameterNames = {'LETd^{ref}','f_{diff}'}; %When optimizing to a reference, one might consider using a quadratic relationship with a non-linear optimizer
        parameterTypes = {'LETd',{'Linear','Quadratic'}};
    end
    
    properties
        parameters = {0,1};        
        penalty = 1;
    end
    
    methods 
        function obj = matRad_MeanLETd(penalty,LETdMeanRef,fDiff)
           
            % if we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@LETdObjectives.matRad_LETdObjective(inputStruct);
            
            if ~initFromStruct
                if nargin < 3 || ~ischar(fDiff)
                    fDiff = 'Linear';
                end
                
                fDiffIx = find(strcmp(fDiff,obj.parameterTypes{2}));
                
                if isempty(fDiffIx) || numel(fDiffIx) > 1
                    fDiffIx = 1;
                    matRad_cfg = MatRad_Config.instance();                    
                    matRad_cfg.dispWarning('Mean LET difference function can only be %s! Using %s difference.', strjoin(obj.parameterTypes{2},' or '), obj.parameterTypes{2}{fDiffIx});
                end
                
                obj.parameters{2} = fDiffIx;


                if nargin >= 2 && isscalar(LETdMeanRef)
                    obj.parameters{1} = LETdMeanRef;
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
        function fLETd = computeLETdObjectiveFunction(obj,LETd)
            switch obj.parameters{2}
                case 1
                    fLETd =  obj.objectiveLinearDiff(LETd);
                case 2
                    fLETd =  obj.objectiveQuadraticDiff(LETd);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid setting for %s in Mean LET Objective!',obj.parameterNames{2});  
            end
        end
        
        %% Calculates the Objective Function gradient
        function fLETdGrad   = computeLETdObjectiveGradient(obj,LETd)
            switch obj.parameters{2}
                case 1
                    fLETdGrad = obj.gradientLinearDiff(LETd);
                case 2
                    fLETdGrad = obj.gradientQuadraticDiff(LETd);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid setting for %s in Mean Dose Objective!',obj.parameterNames{2});  
            end
        end
    end

    methods (Access = protected)
        function fLETd = objectiveQuadraticDiff(obj,LETd)
            fLETd = (mean(LETd(:)) - obj.parameters{1})^2;
        end

        function fLETdGrad = gradientQuadraticDiff(obj,LETd)
            fLETdGrad = 2*(mean(LETd(:))-obj.parameters{1}) * ones(size(LETd(:)))/numel(LETd);
        end

        function fLETd = objectiveLinearDiff(obj,LETd)
            fLETd = abs(mean(LETd(:)) - obj.parameters{1});
        end

        function fLETdGrad = gradientLinearDiff(obj,LETd)
            fLETdGrad = (1/numel(LETd))*sign(LETd(:)-obj.parameters{1});
        end
    end
    
end

