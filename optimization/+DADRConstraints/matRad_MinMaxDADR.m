classdef matRad_MinMaxDADR < DADRConstraints.matRad_DADRConstraint
    % matRad_MinMaxDADR Implements a MinMaxDADR constraint
    %   See matRad_DADRConstraint for interface description
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
        name = 'Min/Max DADR constraint';
        parameterNames = {'DADR^{min}', 'DADR^{max}','method'};
        parameterTypes = {'DADR','DADR',{'approx','voxelwise'}};
    end
    
    properties
        parameters = {40,Inf,1}; %DADR >= 40 Gy/s may be a good example value for FLASH
        epsilon = 1e-3; %slack parameter for the logistic approximation
    end
    
    methods
        function obj = matRad_MinMaxDADR(minDADR,maxDADR,method)
            
            %If we have a struct in first argument
            if nargin == 1 && isstruct(minDADR)
                inputStruct = minDADR;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DADRConstraints.matRad_DADRConstraint(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                
                if nargin < 3 || ~ischar(method)
                    method = 'approx';
                end
                
                methodIx = find(strcmp(method,obj.parameterTypes{3}));
                
                if isempty(methodIx) || numel(methodIx) > 1
                    methodIx = 1;
                    msg = ['Dose Constraint method can only be ', strjoin(obj.parameterTypes{3},' or '), '! Using method ''', obj.parameterTypes{3}{methodIx}, '''.'];
                    warning(msg);
                end
                
                obj.parameters{3} = methodIx;
                
                if nargin >= 2 && isscalar(maxDADR)
                    obj.parameters{2} = maxDADR;
                end
                
                if nargin >= 1 && isscalar(minDADR)
                    obj.parameters{1} = minDADR;
                end
            end
        end
        
        %Overloads the struct function to add constraint specific
        %parameters
        function s = struct(obj)
            s = struct@DADRConstraints.matRad_DADRConstraint(obj);
            s.epsilon = obj.epsilon;
        end
        
        function cu = upperBounds(obj,n)
            switch obj.parameters{3}
                case 1 %logsumexp approx
                    %Validate parameters
                    if obj.parameters{1} <= 0 && isinf(obj.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                        cu = [];
                    elseif obj.parameters{2} == Inf %Only min dose
                        cu = Inf;
                    elseif obj.parameters{1} <= 0 %Only max dose
                        cu = obj.parameters{2};
                    else %both are set sensible
                        cu = [Inf; obj.parameters{2}];
                    end
                case 2 %voxelwise
                    cu = obj.parameters{2}*ones(n,1);
                otherwise
                    error(['Min/max dose constraint evaluation method not known!']);
            end
            %cu = [Inf; obj.parameters{2}];
        end
        function cl = lowerBounds(obj,n)
            switch obj.parameters{3}
                case 1 %logsumexp approx
                    if obj.parameters{1} <= 0 && isinf(obj.parameters{2})
                        cl = [];
                    elseif obj.parameters{2} == Inf
                        cl = obj.parameters{1};
                    elseif obj.parameters{1} <= 0
                        cl = 0;
                    else
                        cl = [obj.parameters{1}; 0];
                    end
                case 2
                    cl = obj.parameters{1}*ones(n,1);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Min/max dose constraint evaluation method not known!');
            end
        end
        
        function jStruct = getDADRConstraintJacobianStructure(obj,n)
            switch obj.parameters{3}
                case 1 %logsumexp approx
                    %Validate parameters
                    if obj.parameters{1} <= 0 && isinf(obj.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                        jStruct = ones(n,0);
                    elseif obj.parameters{1} > 0 && isfinite(obj.parameters{2}) %both are set sensible
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
        function cDose = computeDADRConstraintFunction(obj,DADR)
            switch obj.parameters{3}
                case 1 %logsumexp approx
                    cDose = obj.computeDADRConstraintFunctionLogSumExp(DADR);
                case 2
                    cDose = obj.computeDADRConstraintFunctionVoxelwise(DADR);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Min/max dose constraint evaluation method not known!');
            end
        end
        
        %% Calculates the Constraint jacobian
        function cDoseJacob  = computeDADRConstraintJacobian(obj,DADR)
            switch obj.parameters{3}
                case 1 %logsumexp approx
                    cDoseJacob = obj.computeDADRConstraintJacobianLogSumExp(DADR);
                case 2
                    cDoseJacob = obj.computeDADRConstraintJacobianVoxelwise(DADR);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Min/max dose constraint evaluation method not known!');
            end
        end
    end
    
    methods (Access = private)
        % LogSumExp Approximation
        function cDose = computeDADRConstraintFunctionLogSumExp(obj,DADR)
            dose_min = min(DADR);
            dose_max = max(DADR);
            
            %Validate parameters
            if obj.parameters{1} <= 0 && isinf(obj.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                cDose = [];
            elseif obj.parameters{2} == Inf %Only min dose
                cDose = dose_min - obj.epsilon * log( sum(exp((dose_min - DADR)/obj.epsilon)));
            elseif obj.parameters{1} <= 0 %Only max dose
                cDose = dose_max + obj.epsilon * log( sum(exp((DADR - dose_max)/obj.epsilon)));
            else %both are set sensible
                cDose(2,1) = dose_max + obj.epsilon * log( sum(exp((DADR - dose_max)/obj.epsilon)));
                cDose(1,1) = dose_min - obj.epsilon * log( sum(exp((dose_min - DADR)/obj.epsilon)));
            end
            
        end
        function cDoseJacob  = computeDADRConstraintJacobianLogSumExp(obj,DADR)
            %Validate parameters
            if obj.parameters{1} <= 0 && isinf(obj.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                cDoseJacob = [];
            elseif obj.parameters{2} == Inf %Only min dose
                cDoseJacob(:,1) = exp( (min(DADR)-DADR)/obj.epsilon );
                cDoseJacob(:,1) = cDoseJacob(:,1)/sum(cDoseJacob(:,1));
            elseif obj.parameters{1} <= 0 %Only max dose
                cDoseJacob(:,1) = exp( (DADR-max(DADR))/obj.epsilon );
                cDoseJacob(:,1) = cDoseJacob(:,1)/sum(cDoseJacob(:,1));
            else %both are set sensible
                cDoseJacob(:,1) = exp( (min(DADR)-DADR)/obj.epsilon );
                cDoseJacob(:,1) = cDoseJacob(:,1)/sum(cDoseJacob(:,1));
                
                cDoseJacob(:,2) = exp( (DADR-max(DADR))/obj.epsilon );
                cDoseJacob(:,2) = cDoseJacob(:,2)/sum(cDoseJacob(:,2));
            end
            
            
        end
        
        %Exact voxel-wise
        function cDose = computeDADRConstraintFunctionVoxelwise(obj,DADR)
            cDose = DADR;
        end
        function cDoseJacob  = computeDADRConstraintJacobianVoxelwise(obj,DADR)
            cDoseJacob = speye(numel(DADR),numel(DADR));
        end
    end
    
end


