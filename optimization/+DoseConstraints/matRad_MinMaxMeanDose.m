classdef matRad_MinMaxMeanDose < DoseConstraints.matRad_DoseConstraint
    % matRad_MinMaxMeanDose Implements a MinMaxMeanDose constraint
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
        name = 'mean dose constraint';
        parameterNames = {'\mu_d^{min}', '\mu_d^{max}'};
        %parameterIsDose = logical([1 1]);
        parameterTypes = {'dose','dose'};
    end
    
    properties
        parameters = {0,30};
    end
    
    methods
        function this = matRad_MinMaxMeanDose(minMeanDose,maxMeanDose)
            
            %If we have a struct in first argument
            if nargin == 1 && isstruct(minMeanDose)
                initFromStruct = true;
                inputStruct = minMeanDose;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            % Call Superclass Constructor (for struct initialization)
            this@DoseConstraints.matRad_DoseConstraint(inputStruct);
            
            % now handle initialization from other parameters
            if ~initFromStruct
                
                if nargin == 2 && isscalar(maxMeanDose)
                    this.parameters{2} = maxMeanDose;
                end
                
                if nargin >= 1 && isscalar(minMeanDose)
                    this.parameters{1} = minMeanDose;
                end
                
            end
        end
        
        function cu = upperBounds(this,n)
            cu = this.parameters{2};
        end
        function cl = lowerBounds(this,n)
            cl = this.parameters{1};
        end
        
        %Overloads the struct function to add constraint specific
        %parameters
        function s = struct(this)
            s = struct@DoseConstraints.matRad_DoseConstraint(this);
            %Nothing to do here...
        end
        
        
        %% Calculates the Constraint Function value
        function cDose = computeDoseConstraintFunction(this,dose)
            cDose = mean(dose);
        end
        
        %% Calculates the Constraint jacobian
        function cDoseJacob  = computeDoseConstraintJacobian(this,dose)
            cDoseJacob = ones(numel(dose),1)./numel(dose);
        end
    end
    
end


