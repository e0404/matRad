classdef matRad_MedianDose < DoseObjectives.matRad_DoseObjective
% matRad_SquaredMeanDose implements a penalized squared mean dose objective
%   See matRad_DoseObjective for interface description
%
% References
%   matRad_MeanDose
%   matRad_SquaredDeviation
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
        name = 'Median Dose';
        parameterNames = {'d^{ref}','f_{diff}','n'}; 
        parameterTypes = {'dose','Linear','numeric'};
    end
    
    properties
        parameters = {0,1,0};        
        penalty = 1;
    end
    
    methods 
        function obj = matRad_MedianDose(penalty,dMeanRef,n,fDiff)
           
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
                if nargin < 4 || ~ischar(fDiff)
                    fDiff = 'Linear';
                end
                
                fDiffIx = find(strcmp(fDiff,obj.parameterTypes{3}));
                
                if isempty(fDiffIx) || numel(fDiffIx) > 1
                    fDiffIx = 1;
                    matRad_cfg = MatRad_Config.instance();                    
                    matRad_cfg.dispWarning('Mean dose difference function can only be %s! Using %s difference.', strjoin(obj.parameterTypes{2},' or '), obj.parameterTypes{2}{fDiffIx});
                end
                
                obj.parameters{3} = fDiffIx;
                
                if nargin >= 3 && isscalar(n)
                    obj.parameters{2} = n;
                end

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
        function fDose = computeDoseObjectiveFunction(obj,dose)
            switch obj.parameters{2}
                case 1
                    fDose =  obj.objectiveMedian(dose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid setting for %s in Mean Dose Objective!',obj.parameterNames{2});  
            end
        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad   = computeDoseObjectiveGradient(obj,dose)
            switch obj.parameters{2}
                case 1
                    fDoseGrad = obj.gradientMedian(dose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid setting for %s in Mean Dose Objective!',obj.parameterNames{2});  
            end
        end
    end

     methods (Access = protected)
         function fDose = objectiveMedian(obj,dose)
             obj.parameters{2} = numel(dose(:));
            if mod(obj.parameters{2},2) == 0
                sortDose = sort(dose(:));
                n_2 = sortDose(1,obj.parameters{2}/2);
                n_21 = sortDose(1,(obj.parameters{2}/2+1));
                fDose = 0.5*(n_2 + n_21);

            elseif mod(obj.parameters{2},2) == 1
                sortDose = sort(dose(:));
                fDose = sortDose(1,((obj.parameters{2}+1)/2));
            end
        end

        function fDoseGrad = gradientMedian(obj,dose)
           obj.parameters{2} = numeld(dose(:));
           if mod(obj.parameters{2},2) == 0
               % fDoseGrad = 
           elseif mod(obj.parameters{2},2) == 1
           end
          
        end

    end
    
end


