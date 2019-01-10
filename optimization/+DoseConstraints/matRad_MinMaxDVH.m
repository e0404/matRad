classdef matRad_MinMaxDVH < DoseConstraints.matRad_DoseConstraint
    %MATRAD_DOSEOBJECTIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        name = 'DVH constraint';
        parameterNames = {'d^{ref}', 'V^{min}', 'V^{max}'};
        %parameterIsDose = logical([1 0 0]);
        parameterTypes = {'dose','numeric','numeric'};
    end
    
    properties
        voxelScalingRatio = 1;
        referenceScalingVal = 0.01;
        parameters = {30,0,100};
    end
    
    methods
        function obj = matRad_MinMaxDVH(dRef,vMin,vMax)
            
            %If we have a struct in first argument
            if nargin == 1 && isstruct(dRef)
                inputStruct = dRef;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DoseConstraints.matRad_DoseConstraint(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct                
                if nargin == 3 && isscalar(vMax)
                    obj.parameters{3} = vMax;
                end
                
                if nargin >= 1 && isscalar(dRef)
                    obj.parameters{1} = dRef;
                end
                
                if nargin >= 2 && isscalar(vMin)
                    obj.parameters{2} = vMin;
                end
            end
        end
        
        %Overloads the struct function to add constraint specific
        %parameters
        function s = struct(obj)
            s = struct@DoseConstraints.matRad_DoseConstraint(obj);
            s.voxelScalingRatio = 1;
            s.referenceScalingVal = 0.01;
        end
        
        function cu = upperBounds(obj,n)
            cu = obj.parameters{3} / 100;
        end
        function cl = lowerBounds(obj,n)
            cl = obj.parameters{2} / 100;
        end
        %% Calculates the Constraint Function value
        function cDose = computeDoseConstraintFunction(obj,dose)
            
            %Fast DVH point calculation
            cDose = sum(dose >= obj.parameters{1})/numel(dose);
            
            %cDose = 100 * cDose; %In Percent
            
            % alternative constraint calculation 3/4 %
            % % get reference Volume
            % refVol = cst{j,6}(k).volume/100;
            %
            % % calc deviation
            % deviation = d_i - d_ref;
            %
            % % calc d_ref2: V(d_ref2) = refVol
            % d_ref2 = matRad_calcInversDVH(refVol,d_i);
            %
            % % apply lower and upper dose limits
            % if isequal(cst{j,6}(k).type, 'max DVH constraint')
            %    deviation(d_i < d_ref | d_i > d_ref2) = 0;
            % elseif isequal(cst{j,6}(k).type, 'min DVH constraint')
            %    deviation(d_i > d_ref | d_i < d_ref2) = 0;
            % end
            %
            % %c = sum(deviation);                              % linear deviation
            % %c = deviation'*deviation;                        % square devioation
            % c = (1/size(cst{j,4},1))*(deviation'*deviation); % square deviation with normalization
            % %c = (deviation).^2'*(deviation).^2;               % squared square devioation
            % alternative constraint calculation 3/4 %
        end
        
        %% Calculates the Constraint jacobian
        function cDoseJacob  = computeDoseConstraintJacobian(obj,dose)
            %logistic approximation
            
            %Do we really need to sort two times?
            dose_sort = sort(dose);
            
            % calculate scaling
            NoVoxels     = max(obj.voxelScalingRatio*numel(dose),10);
            absDiffsort  = sort(abs(obj.parameters{1} - dose_sort));
            deltaDoseMax = absDiffsort(ceil(NoVoxels/2));
            
            % calclulate DVHC scaling
            DVHCScaling = min((log(1/obj.referenceScalingVal-1))/(2*deltaDoseMax),250);
            
            d_diff = dose - obj.parameters{1};
            
            cDoseJacob = (2/numel(dose))*DVHCScaling*exp(2*DVHCScaling*d_diff)./(exp(2*DVHCScaling*d_diff)+1).^2;
            
            % alternative constraint calculation 4/4 %
            % % get reference Volume
            % refVol = cst{j,6}(k).volume/100;
            %
            % % calc deviation
            % deviation = d_i - d_ref;
            %
            % % calc d_ref2: V(d_ref2) = refVol
            % d_ref2 = matRad_calcInversDVH(refVol,d_i);
            %
            % % apply lower and upper dose limits
            % if isequal(cst{j,6}(k).type, 'max DVH constraint')
            %      deviation(d_i < d_ref | d_i > d_ref2) = 0;
            % elseif isequal(cst{j,6}(k).type, 'min DVH constraint')
            %      deviation(d_i > d_ref | d_i < d_ref2) = 0;
            % end
            %
            % %jacobVec = ones(size(cst{j,4}));             % linear deviation
            % %jacobVec = 2*deviation;                      % square deviation
            % jacobVec = (1/size(cst{j,4},1))*2*deviation; % square deviation with normalization
            % %jacobVec = 4*(deviation).^3;                  % squared square devioation
            % alternative constraint calculation 4/4 %
        end
    end
    
end


