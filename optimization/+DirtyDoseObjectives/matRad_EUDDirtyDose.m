classdef matRad_EUDDirtyDose < DirtyDoseObjectives.matRad_DirtyDoseObjective
% matRad_EUDDirtyDose Implements a penalized equivalent uniform dirty dose objective
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
        name = 'EUD';
        parameterNames = {'EUD^{ref}', 'k'};
        parameterTypes = {'dirtyDose','numeric'};
    end
    
    properties
        parameters = {0, 3.5};
        penalty = 1;
    end
    
    methods
        function obj = matRad_EUDDirtyDose(penalty,eudRef, eudExponent)
            %If we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DirtyDoseObjectives.matRad_DirtyDoseObjective(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                if nargin >= 3 && isscalar(eudExponent)
                    obj.parameters{2} = eudExponent;
                end
                
                if nargin >= 2 && isscalar(eudRef)
                    obj.parameters{1} = eudRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
        
        %% Calculates the Objective Function value
        function fDirtyDose = computeDirtyDoseObjectiveFunction(obj,dirtyDose)
            % get exponent for EUD
            k = obj.parameters{2};
            
            % calculate power sum
            powersum = sum(dirtyDose.^k);
           
            %Calculate objective
            fDirtyDose = (nthroot(powersum/numel(dirtyDose),k) - obj.parameters{1})^2;
        end
        
        %% Calculates the Objective Function gradient
        function fDirtyDoseGrad  = computeDirtyDoseObjectiveGradient(obj,dirtyDose)
            % get exponent for EUD
            k = obj.parameters{2};
            
            %numerical stability
            dirtyDose(dirtyDose == 0) = 0.001;
            
            % calculate power sum
            powersum = sum(dirtyDose.^k);

            %derivatives = nthroot(1/numel(dose),k) * powersum^((1-k)/k) * (dose.^(k-1));
            fDirtyDoseGrad = 2 * nthroot(1/numel(dirtyDose),k) * powersum^((1-k)/k) * (dirtyDose.^(k-1)) .* (nthroot(powersum/numel(dirtyDose),k) - obj.parameters{1});
            %end
            if any(~isfinite(fDirtyDoseGrad)) % check for inf and nan for numerical stability
                error('EUD computation failed. Reduce exponent to resolve numerical problems.');
            end
        end
    end
    
end

