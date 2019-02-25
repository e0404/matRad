classdef matRad_MinDVH < DoseObjectives.matRad_DoseObjective
    %MATRAD_DOSEOBJECTIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        name = 'Min DVH';
        parameterNames = {'d', 'V^{min}'};
        parameterTypes = {'dose','numeric'};
    end
    
    properties
        parameters = {60,95};
        penalty = 1;
    end
    
    methods
        function obj = matRad_MinDVH(penalty,dRef,vMinPercent)
            
            % if we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DoseObjectives.matRad_DoseObjective(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                if nargin >= 3 && isscalar(vMinPercent)
                    obj.parameters{2} = vMinPercent;
                end
                
                if nargin >= 2 && isscalar(dRef)
                    obj.parameters{1} = dRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
            
        end        
        %% Calculates the Objective Function value
        function fDose = computeDoseObjectiveFunction(obj,dose)                       
            % get reference Volume
            refVol = obj.parameters{2}/100;
            
            % calc deviation
            deviation = dose - obj.parameters{1};

            % calc d_ref2: V(d_ref2) = refVol
            d_ref2 = matRad_calcInversDVH(refVol,dose);

            
            deviation(dose > obj.parameters{1} | dose < d_ref2) = 0;
   
            % claculate objective function
            fDose = (obj.penalty/numel(dose))*(deviation'*deviation);
        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad   = computeDoseObjectiveGradient(obj,dose)
            % get reference Volume
            refVol = obj.parameters{2}/100;
            
            % calc deviation
            deviation = dose - obj.parameters{1};
            
            % calc d_ref2: V(d_ref2) = refVol
            d_ref2 = matRad_calcInversDVH(refVol,dose);
            
            deviation(dose > obj.parameters{1} | dose < d_ref2) = 0;

            % calculate delta
            fDoseGrad = 2 * (obj.penalty/numel(dose))*deviation;
        end
    end
    
end