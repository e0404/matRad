classdef matRad_MinMaxDose < DoseConstraints.matRad_DoseConstraint
    %MATRAD_DOSEOBJECTIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    % use log sum exp approximation, see appendix A in
    % http://scitation.aip.org/content/aapm/journal/medphys/41/8/10.1118/1.4883837
    
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
        function obj = matRad_MinMaxDose(minDose,maxDose,method)
            
            %If we have a struct in first argument
            if nargin == 1 && isstruct(minDose)
                inputStruct = minDose;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DoseConstraints.matRad_DoseConstraint(inputStruct);
            
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
                
                if nargin >= 2 && isscalar(maxDose)
                    obj.parameters{2} = maxDose;
                end
                
                if nargin >= 1 && isscalar(minDose)
                    obj.parameters{1} = minDose;
                end
            end
        end
        
        %Overloads the struct function to add constraint specific
        %parameters
        function s = struct(obj)
            s = struct@DoseConstraints.matRad_DoseConstraint(obj);
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
                    error(['Min/max dose constraint evaluation method ''' obj.method ''' not known!']);
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
                    error(['Min/max dose constraint evaluation method ''' obj.method ''' not known!']);
            end
        end
        
        function jStruct = getDoseConstraintJacobianStructure(obj,n)
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
                    error(['Min/max dose constraint evaluation method ''' obj.method ''' not known!']);
            end
            
        end
        
        %% Calculates the Constraint Function value
        function cDose = computeDoseConstraintFunction(obj,dose)
            %cDose(2) = dose_max + obj.epsilon * log( sum(exp((dose - dose_max)/obj.epsilon)));
            %cDose(1) = dose_min - obj.epsilon * log( sum(exp((dose_min - dose)/obj.epsilon)));
            switch obj.parameters{3}
                case 1 %logsumexp approx
                    cDose = obj.computeDoseConstraintFunctionLogSumExp(dose);
                case 2
                    cDose = obj.computeDoseConstraintFunctionVoxelwise(dose);
                otherwise
                    error(['Min/max dose constraint evaluation method ''' obj.method ''' not known!']);
            end
        end
        
        %% Calculates the Constraint jacobian
        function cDoseJacob  = computeDoseConstraintJacobian(obj,dose)
            switch obj.parameters{3}
                case 1 %logsumexp approx
                    cDoseJacob = obj.computeDoseConstraintJacobianLogSumExp(dose);
                case 2
                    cDoseJacob = obj.computeDoseConstraintJacobianVoxelwise(dose);
                otherwise
                    error(['Min/max dose constraint evaluation method ''' obj.method ''' not known!']);
            end
        end
    end
    
    methods (Access = private)
        % LogSumExp Approximation
        function cDose = computeDoseConstraintFunctionLogSumExp(obj,dose)
            dose_min = min(dose);
            dose_max = max(dose);
            
            %Validate parameters
            if obj.parameters{1} <= 0 && isinf(obj.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                cDose = [];
            elseif obj.parameters{2} == Inf %Only min dose
                cDose = dose_min - obj.epsilon * log( sum(exp((dose_min - dose)/obj.epsilon)));
            elseif obj.parameters{1} <= 0 %Only max dose
                cDose = dose_max + obj.epsilon * log( sum(exp((dose - dose_max)/obj.epsilon)));
            else %both are set sensible
                cDose(2) = dose_max + obj.epsilon * log( sum(exp((dose - dose_max)/obj.epsilon)));
                cDose(1) = dose_min - obj.epsilon * log( sum(exp((dose_min - dose)/obj.epsilon)));
            end
            
        end
        function cDoseJacob  = computeDoseConstraintJacobianLogSumExp(obj,dose)
            %Validate parameters
            if obj.parameters{1} <= 0 && isinf(obj.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                cDoseJacob = [];
            elseif obj.parameters{2} == Inf %Only min dose
                cDoseJacob(:,1) = exp( (min(dose)-dose)/obj.epsilon );
                cDoseJacob(:,1) = cDoseJacob(:,1)/sum(cDoseJacob(:,1));
            elseif obj.parameters{1} <= 0 %Only max dose
                cDoseJacob(:,1) = exp( (dose-max(dose))/obj.epsilon );
                cDoseJacob(:,1) = cDoseJacob(:,1)/sum(cDoseJacob(:,1));
            else %both are set sensible
                cDoseJacob(:,1) = exp( (min(dose)-dose)/obj.epsilon );
                cDoseJacob(:,1) = cDoseJacob(:,1)/sum(cDoseJacob(:,1));
                
                cDoseJacob(:,2) = exp( (dose-max(dose))/obj.epsilon );
                cDoseJacob(:,2) = cDoseJacob(:,2)/sum(cDoseJacob(:,2));
            end
            
            
        end
        
        %Exact voxel-wise
        function cDose = computeDoseConstraintFunctionVoxelwise(obj,dose)
            cDose = dose;
        end
        function cDoseJacob  = computeDoseConstraintJacobianVoxelwise(obj,dose)
            cDoseJacob = speye(numel(dose),numel(dose));
        end
    end
    
end


