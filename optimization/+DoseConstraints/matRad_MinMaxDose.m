classdef matRad_MinMaxDose < DoseConstraints.matRad_DoseConstraint
    % matRad_MinMaxDose Implements a MinMaxDose constraint
    %   See matRad_DoseConstraint for interface description
    %
    % use log sum exp approximation, see appendix A in
    % http://scitation.aip.org/content/aapm/journal/medphys/41/8/10.1118/1.4883837
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
        name = 'Min/Max dose constraint';
        parameterNames = {'d^{min}', 'd^{max}','method'};
        parameterTypes = {'dose','dose',{'approx','voxelwise'}};
    end
    
    properties
        parameters = {0,30,1};
        epsilon = 1e-3; %slack parameter for the logistic approximation
    end
    
    methods
        function this = matRad_MinMaxDose(minDose,maxDose,method)
            
            %If we have a struct in first argument
            if nargin == 1 && isstruct(minDose)
                inputStruct = minDose;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            this@DoseConstraints.matRad_DoseConstraint(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                
                if nargin < 3 || ~ischar(method)
                    method = 'approx';
                end
                
                methodIx = find(strcmp(method,this.parameterTypes{3}));
                
                if isempty(methodIx) || numel(methodIx) > 1
                    methodIx = 1;
                    msg = ['Dose Constraint method can only be ', strjoin(this.parameterTypes{3},' or '), '! Using method ''', this.parameterTypes{3}{methodIx}, '''.'];
                    warning(msg);
                end
                
                this.parameters{3} = methodIx;
                
                if nargin >= 2 && isscalar(maxDose)
                    this.parameters{2} = maxDose;
                end
                
                if nargin >= 1 && isscalar(minDose)
                    this.parameters{1} = minDose;
                end
            end
        end
        
        %Overloads the struct function to add constraint specific
        %parameters
        function s = struct(this)
            s = struct@DoseConstraints.matRad_DoseConstraint(this);
            s.epsilon = this.epsilon;
        end
        
        function cu = upperBounds(this,n)
            switch this.parameters{3}
                case 1 %logsumexp approx
                    %Validate parameters
                    if this.parameters{1} <= 0 && isinf(this.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                        cu = [];
                    elseif this.parameters{2} == Inf %Only min dose
                        cu = Inf;
                    elseif this.parameters{1} <= 0 %Only max dose
                        cu = this.parameters{2};
                    else %both are set sensible
                        cu = [Inf; this.parameters{2}];
                    end
                case 2 %voxelwise
                    cu = this.parameters{2}*ones(n,1);
                otherwise
                    error(['Min/max dose constraint evaluation method not known!']);
            end
            %cu = [Inf; this.parameters{2}];
        end
        function cl = lowerBounds(this,n)
            switch this.parameters{3}
                case 1 %logsumexp approx
                    if this.parameters{1} <= 0 && isinf(this.parameters{2})
                        cl = [];
                    elseif this.parameters{2} == Inf
                        cl = this.parameters{1};
                    elseif this.parameters{1} <= 0
                        cl = 0;
                    else
                        cl = [this.parameters{1}; 0];
                    end
                case 2
                    cl = this.parameters{1}*ones(n,1);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Min/max dose constraint evaluation method not known!');
            end
        end
        
        function jStruct = getDoseConstraintJacobianStructure(this,n)
            switch this.parameters{3}
                case 1 %logsumexp approx
                    %Validate parameters
                    if this.parameters{1} <= 0 && isinf(this.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                        jStruct = ones(n,0);
                    elseif this.parameters{1} > 0 && isfinite(this.parameters{2}) %both are set sensible
                        jStruct = ones(n,2);
                    else %Only min or max dose
                        jStruct = ones(n,1);
                    end
                    %jStruct = ones(n,2);
                case 2
                    jStruct = speye(n);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Min/max dose constraint evaluation method not known!');
            end
            
        end
        
        %% Calculates the Constraint Function value
        function cDose = computeDoseConstraintFunction(this,dose)
            %cDose(2) = dose_max + this.epsilon * log( sum(exp((dose - dose_max)/this.epsilon)));
            %cDose(1) = dose_min - this.epsilon * log( sum(exp((dose_min - dose)/this.epsilon)));
            switch this.parameters{3}
                case 1 %logsumexp approx
                    cDose = this.computeDoseConstraintFunctionLogSumExp(dose);
                case 2
                    cDose = this.computeDoseConstraintFunctionVoxelwise(dose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Min/max dose constraint evaluation method not known!');
            end
        end
        
        %% Calculates the Constraint jacobian
        function cDoseJacob  = computeDoseConstraintJacobian(this,dose)
            switch this.parameters{3}
                case 1 %logsumexp approx
                    cDoseJacob = this.computeDoseConstraintJacobianLogSumExp(dose);
                case 2
                    cDoseJacob = this.computeDoseConstraintJacobianVoxelwise(dose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Min/max dose constraint evaluation method not known!');
            end
        end
    end
    
    methods (Access = private)
        % LogSumExp Approximation
        function cDose = computeDoseConstraintFunctionLogSumExp(this,dose)
            dose_min = min(dose);
            dose_max = max(dose);
            
            %Validate parameters
            if this.parameters{1} <= 0 && isinf(this.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                cDose = [];
            elseif this.parameters{2} == Inf %Only min dose
                cDose = dose_min - this.epsilon * log( sum(exp((dose_min - dose)/this.epsilon)));
            elseif this.parameters{1} <= 0 %Only max dose
                cDose = dose_max + this.epsilon * log( sum(exp((dose - dose_max)/this.epsilon)));
            else %both are set sensible
                cDose(2,1) = dose_max + this.epsilon * log( sum(exp((dose - dose_max)/this.epsilon)));
                cDose(1,1) = dose_min - this.epsilon * log( sum(exp((dose_min - dose)/this.epsilon)));
            end
            
        end
        function cDoseJacob  = computeDoseConstraintJacobianLogSumExp(this,dose)
            %Validate parameters
            if this.parameters{1} <= 0 && isinf(this.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                cDoseJacob = [];
            elseif this.parameters{2} == Inf %Only min dose
                cDoseJacob(:,1) = exp( (min(dose)-dose)/this.epsilon );
                cDoseJacob(:,1) = cDoseJacob(:,1)/sum(cDoseJacob(:,1));
            elseif this.parameters{1} <= 0 %Only max dose
                cDoseJacob(:,1) = exp( (dose-max(dose))/this.epsilon );
                cDoseJacob(:,1) = cDoseJacob(:,1)/sum(cDoseJacob(:,1));
            else %both are set sensible
                cDoseJacob(:,1) = exp( (min(dose)-dose)/this.epsilon );
                cDoseJacob(:,1) = cDoseJacob(:,1)/sum(cDoseJacob(:,1));
                
                cDoseJacob(:,2) = exp( (dose-max(dose))/this.epsilon );
                cDoseJacob(:,2) = cDoseJacob(:,2)/sum(cDoseJacob(:,2));
            end
            
            
        end
        
        %Exact voxel-wise
        function cDose = computeDoseConstraintFunctionVoxelwise(this,dose)
            cDose = dose;
        end
        function cDoseJacob  = computeDoseConstraintJacobianVoxelwise(this,dose)
            cDoseJacob = speye(numel(dose),numel(dose));
        end
    end
    
end


