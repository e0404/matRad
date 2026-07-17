 classdef matRad_MinMaxDVH < DoseConstraints.matRad_DoseConstraint
    % matRad_MinMaxDVH Implements a MinMaxDVH constraint
    %   See matRad_DoseConstraint for interface description
    %
    % The constraint enforces a DVH point V(d^{ref}) - the volume fraction of
    % the structure receiving at least the reference dose d^{ref} - to lie
    % between V^{min} and V^{max}.
    %
    % The constraint value is the exact DVH point (a Heaviside/step count). Its
    % exact derivative is a Dirac delta (zero everywhere except at the
    % threshold) and is therefore useless for gradient-based optimization. For
    % the Jacobian we hence approximate the Heaviside step by a logistic
    % (sigmoid) function whose steepness is re-derived in every iteration - the
    % dose distribution within a structure changes during optimization, so the
    % steepness is adapted such that a "meaningful" gradient is provided for the
    % voxels close to the threshold (those that can realistically be pushed
    % across the DVH point). The steepness is controlled by two parameters:
    %
    %   voxelScalingRatio   - WIDTH of the sensitivity band around the
    %                         threshold, expressed via a voxel count. The band
    %                         half-width deltaDoseMax is the dose distance
    %                         enclosing the closest voxelScalingRatio*N/2 voxels
    %                         (default 1 -> median distance, i.e. the closest
    %                         half of all voxels). Smaller values narrow the
    %                         band (sharper sigmoid, fewer contributing voxels)
    %                         and can be useful for OARs, where only few voxels
    %                         violate the constraint.
    %   referenceScalingVal - HOW step-like the sigmoid is at the edge of that
    %                         band: at +/- deltaDoseMax the sigmoid must be
    %                         within referenceScalingVal of its asymptote
    %                         (default 0.01 -> saturated to within 1%). Smaller
    %                         values make the transition crisper. Largely
    %                         empirical; it can make sense to reduce it as the
    %                         optimization settles close to a constraint value.
    %
    % See the Fluence Optimization section of the documentation for details.
    %
    % References
    %   -
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2020-2026 the matRad development team.
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Constant)
        name = 'DVH constraint';
        parameterNames = {'d^{ref}', 'V^{min}', 'V^{max}'};
        %parameterIsDose = logical([1 0 0]);
        parameterTypes = {'dose','numeric','numeric'};
    end
    
    properties
        voxelScalingRatio = 1;      %width of the logistic sensitivity band (fraction of voxels); see class help
        referenceScalingVal = 0.01; %sigmoid saturation residual at the band edge; see class help
        parameters = {30,0,100};
    end
    
    methods
        function this = matRad_MinMaxDVH(dRef,vMin,vMax)
            
            %If we have a struct in first argument
            if nargin == 1 && isstruct(dRef)
                inputStruct = dRef;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            this@DoseConstraints.matRad_DoseConstraint(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct                
                if nargin == 3 && isscalar(vMax)
                    this.parameters{3} = vMax;
                end
                
                if nargin >= 1 && isscalar(dRef)
                    this.parameters{1} = dRef;
                end
                
                if nargin >= 2 && isscalar(vMin)
                    this.parameters{2} = vMin;
                end
            end
        end
        
        %Overloads the struct function to add constraint specific
        %parameters
        function s = struct(this)
            s = struct@DoseConstraints.matRad_DoseConstraint(this);
            s.voxelScalingRatio = 1;
            s.referenceScalingVal = 0.01;
        end
        
        function cu = upperBounds(this,n)
            cu = this.parameters{3} / 100;
        end
        function cl = lowerBounds(this,n)
            cl = this.parameters{2} / 100;
        end
        %% Calculates the Constraint Function value
        function cDose = computeDoseConstraintFunction(this,dose)
            
            %Fast DVH point calculation
            cDose = sum(dose >= this.parameters{1})/numel(dose);
            
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
        function cDoseJacob  = computeDoseConstraintJacobian(this,dose)
            % Jacobian of the logistic (sigmoid) approximation to the
            % Heaviside DVH count. The sigmoid steepness DVHCScaling is adapted
            % to the current dose distribution every call (see class help).

            %Do we really need to sort two times?
            dose_sort = sort(dose);

            % Target voxel count defining the sensitivity band width. With the
            % default voxelScalingRatio = 1 this is all voxels (floored at 10).
            NoVoxels     = max(this.voxelScalingRatio*numel(dose),10);

            % Voxel doses sorted by absolute distance from the reference dose.
            absDiffsort  = sort(abs(this.parameters{1} - dose_sort));

            % Band half-width: distance enclosing the closest NoVoxels/2 voxels
            % (the /2 makes the symmetric band +/- deltaDoseMax hold ~NoVoxels/2
            % voxels). With the default ratio this is the median distance.
            deltaDoseMax = absDiffsort(min(ceil(NoVoxels/2),numel(dose)));

            % Steepness k so that the sigmoid is within referenceScalingVal of
            % its asymptote at +/- deltaDoseMax, i.e. sigma(deltaDoseMax) =
            % 1 - referenceScalingVal. Capped at 250 to avoid overflow / an
            % overly localized gradient when deltaDoseMax is tiny.
            DVHCScaling = min((log(1/this.referenceScalingVal-1))/(2*deltaDoseMax),250);

            d_diff = dose - this.parameters{1};

            % Per-voxel derivative of the sigmoid: a smooth bump centered at
            % the threshold (2/N * k * sigma * (1 - sigma)).
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


