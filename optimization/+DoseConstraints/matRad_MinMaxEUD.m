classdef matRad_MinMaxEUD < DoseConstraints.matRad_DoseConstraint
    % matRad_MinMaxEUD Implements a MinMaxEUD constraint
    %   See matRad_DoseConstraint for interface description
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
        name = 'EUD constraint';
        parameterNames = {'k','EUD^{min}', 'EUD^{max}'};
        %parameterIsDose = logical([0 1 1]);
        parameterTypes = {'numeric','dose','dose'};
    end
    
    properties
        parameters = {5,0,30};
    end
    
    methods
        function obj = matRad_MinMaxEUD(exponent,eudMin,eudMax)
            %If we have a struct in first argument
            if nargin == 1 && isstruct(exponent)
                inputStruct = exponent;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DoseConstraints.matRad_DoseConstraint(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                
                if nargin == 3 && isscalar(eudMax)
                    obj.parameters{3} = eudMax;
                end
                
                if nargin >= 1 && isscalar(exponent)
                    obj.parameters{1} = exponent;
                end
                
                if nargin >= 2 && isscalar(eudMin)
                    obj.parameters{2} = eudMin;
                end
            end
        end
        
        %Overloads the struct function to add constraint specific
        %parameters
        function s = struct(obj)
            s = struct@DoseConstraints.matRad_DoseConstraint(obj);
            %Nothing to do here...
        end
        
        function cu = upperBounds(obj,n)
            cu = obj.parameters{3};
        end
        function cl = lowerBounds(obj,n)
            cl = obj.parameters{2};
        end
        %% Calculates the Constraint Function value
        function cDose = computeDoseConstraintFunction(obj,dose)
            k = obj.parameters{1};
            cDose = mean(dose.^k)^(1/k);
        end
        
        %% Calculates the Constraint jacobian
        function cDoseJacob  = computeDoseConstraintJacobian(obj,dose)
            k = obj.parameters{1};
            cDoseJacob = nthroot(1/numel(dose),k) * sum(dose.^k)^((1-k)/k) * (dose.^(k-1));
        end
    end
    
end


