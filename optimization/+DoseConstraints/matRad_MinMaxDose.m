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
        function constr = matRad_MinMaxDose(minDose,maxDose,method)
            
            %If we have a struct in first argument
            if nargin == 1 && isstruct(minDose)
                inputStruct = minDose;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            constr@DoseConstraints.matRad_DoseConstraint(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                
                if nargin < 3 || ~ischar(method)
                    method = 'approx';
                end
                
                methodIx = find(strcmp(method,constr.parameterTypes{3}));
                
                if isempty(methodIx) || numel(methodIx) > 1
                    methodIx = 1;
                    msg = ['Dose Constraint method can only be ', strjoin(constr.parameterTypes{3},' or '), '! Using method ''', constr.parameterTypes{3}{methodIx}, '''.'];
                    warning(msg);
                end
                
                constr.parameters{3} = methodIx;
                
                if nargin >= 2 && isscalar(maxDose)
                    constr.parameters{2} = maxDose;
                end
                
                if nargin >= 1 && isscalar(minDose)
                    constr.parameters{1} = minDose;
                end
            end
        end
        
        %Overloads the struct function to add constraint specific
        %parameters
        function s = struct(constr)
            s = struct@DoseConstraints.matRad_DoseConstraint(constr);
            s.epsilon = constr.epsilon;
        end
        
        function cu = upperBounds(constr,n)
            switch constr.parameters{3}
                case 1 %logsumexp approx
                    %Validate parameters
                    if constr.parameters{1} <= 0 && isinf(constr.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                        cu = [];
                    elseif constr.parameters{2} == Inf %Only min dose
                        cu = Inf;
                    elseif constr.parameters{1} <= 0 %Only max dose
                        cu = constr.parameters{2};
                    else %both are set sensible
                        cu = [Inf; constr.parameters{2}];
                    end
                case 2 %voxelwise
                    cu = constr.parameters{2}*ones(n,1);
                otherwise
                    error(['Min/max dose constraint evaluation method not known!']);
            end
            %cu = [Inf; constr.parameters{2}];
        end
        function cl = lowerBounds(constr,n)
            switch constr.parameters{3}
                case 1 %logsumexp approx
                    if constr.parameters{1} <= 0 && isinf(constr.parameters{2})
                        cl = [];
                    elseif constr.parameters{2} == Inf
                        cl = constr.parameters{1};
                    elseif constr.parameters{1} <= 0
                        cl = 0;
                    else
                        cl = [constr.parameters{1}; 0];
                    end
                case 2
                    cl = constr.parameters{1}*ones(n,1);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Min/max dose constraint evaluation method not known!');
            end
        end
        
        function jStruct = getDoseConstraintJacobianStructure(constr,n)
            switch constr.parameters{3}
                case 1 %logsumexp approx
                    %Validate parameters
                    if constr.parameters{1} <= 0 && isinf(constr.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                        jStruct = ones(n,0);
                    elseif constr.parameters{1} > 0 && isfinite(constr.parameters{2}) %both are set sensible
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
        function cDose = computeDoseConstraintFunction(constr,dose)
            %cDose(2) = dose_max + constr.epsilon * log( sum(exp((dose - dose_max)/constr.epsilon)));
            %cDose(1) = dose_min - constr.epsilon * log( sum(exp((dose_min - dose)/constr.epsilon)));
            switch constr.parameters{3}
                case 1 %logsumexp approx
                    cDose = constr.computeDoseConstraintFunctionLogSumExp(dose);
                case 2
                    cDose = constr.computeDoseConstraintFunctionVoxelwise(dose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Min/max dose constraint evaluation method not known!');
            end
        end
        
        %% Calculates the Constraint jacobian
        function cDoseJacob  = computeDoseConstraintJacobian(constr,dose)
            switch constr.parameters{3}
                case 1 %logsumexp approx
                    cDoseJacob = constr.computeDoseConstraintJacobianLogSumExp(dose);
                case 2
                    cDoseJacob = constr.computeDoseConstraintJacobianVoxelwise(dose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Min/max dose constraint evaluation method not known!');
            end
        end
    end
    
    methods (Access = private)
        % LogSumExp Approximation
        function cDose = computeDoseConstraintFunctionLogSumExp(constr,dose)
            dose_min = min(dose);
            dose_max = max(dose);
            
            %Validate parameters
            if constr.parameters{1} <= 0 && isinf(constr.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                cDose = [];
            elseif constr.parameters{2} == Inf %Only min dose
                cDose = dose_min - constr.epsilon * log( sum(exp((dose_min - dose)/constr.epsilon)));
            elseif constr.parameters{1} <= 0 %Only max dose
                cDose = dose_max + constr.epsilon * log( sum(exp((dose - dose_max)/constr.epsilon)));
            else %both are set sensible
                cDose(2,1) = dose_max + constr.epsilon * log( sum(exp((dose - dose_max)/constr.epsilon)));
                cDose(1,1) = dose_min - constr.epsilon * log( sum(exp((dose_min - dose)/constr.epsilon)));
            end
            
        end
        function cDoseJacob  = computeDoseConstraintJacobianLogSumExp(constr,dose)
            %Validate parameters
            if constr.parameters{1} <= 0 && isinf(constr.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                cDoseJacob = [];
            elseif constr.parameters{2} == Inf %Only min dose
                cDoseJacob(:,1) = exp( (min(dose)-dose)/constr.epsilon );
                cDoseJacob(:,1) = cDoseJacob(:,1)/sum(cDoseJacob(:,1));
            elseif constr.parameters{1} <= 0 %Only max dose
                cDoseJacob(:,1) = exp( (dose-max(dose))/constr.epsilon );
                cDoseJacob(:,1) = cDoseJacob(:,1)/sum(cDoseJacob(:,1));
            else %both are set sensible
                cDoseJacob(:,1) = exp( (min(dose)-dose)/constr.epsilon );
                cDoseJacob(:,1) = cDoseJacob(:,1)/sum(cDoseJacob(:,1));
                
                cDoseJacob(:,2) = exp( (dose-max(dose))/constr.epsilon );
                cDoseJacob(:,2) = cDoseJacob(:,2)/sum(cDoseJacob(:,2));
            end
            
            
        end
        
        %Exact voxel-wise
        function cDose = computeDoseConstraintFunctionVoxelwise(constr,dose)
            cDose = dose;
        end
        function cDoseJacob  = computeDoseConstraintJacobianVoxelwise(constr,dose)
            cDoseJacob = speye(numel(dose),numel(dose));
        end
    end
    
end


