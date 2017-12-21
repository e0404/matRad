classdef matRad_MinMaxDVH < DoseConstraints.matRad_DoseConstraint
    %MATRAD_DOSEOBJECTIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        name = 'Min/Max DVH'
    end
    
    properties
        parameters = {'Reference Dose (Max)', 'Reference Min Volume [%]', 'Reference Max Volume [%]'; 60,0,1}
        voxelScalingRatio = 1;
        referenceScalingVal = 0.01;
    end
    
    methods
        function cu = upperBound(obj)
            cu = obj.parameters{2,3};
        end
        function cl = lowerBound(obj)
            cl = obj.parameters{2,2};
        end
        %% Calculates the Objective Function value
        function cDose = computeDoseConstraintFunction(obj,dose)
            cDose = sum(dose >= obj.parameters{2,1})/numel(dose);
            
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
        
        %% Calculates the Objective Function gradient
        function cDoseJacob  = computeDoseConstraintJacobian(obj,dose)
            %logistic approximation
            
            %Do we really need to sort two times?
            dose_sort = sort(dose);
            
            % calculate scaling
            NoVoxels     = max(obj.voxelScalingRatio*numel(dose),10);
            absDiffsort  = sort(abs(obj.parameters{2,1} - dose_sort));
            deltaDoseMax = absDiffsort(ceil(NoVoxels/2));
            
            % calclulate DVHC scaling
            DVHCScaling = min((log(1/obj.referenceScalingVal-1))/(2*deltaDoseMax),250);
                       
            d_diff = dose - obj.parameters{2,1};   
            
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


