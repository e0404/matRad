classdef (Abstract) matRad_DoseObjective
    %MATRAD_DOSEOBJECTIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract, Constant)        
        name
        parameterNames
        parameterIsDose
    end
    
    properties (Constant)
        type = 'Objective';
    end
    
    properties (Abstract, Access = public)
        parameters        
        penalty
    end
    
    methods (Abstract)
        fDose       = computeDoseObjectiveFunction(obj,dose)
        fDoseGrad   = computeDoseObjectiveGradient(obj,dose)
    end        
    
    %Helper methods
    methods (Access = public)
        %Get only the parameters describing some kind of reference dose as
        %numeric array
        function doseParams = getDoseParameters(obj)
            doseParams = [obj.parameters{obj.parameterIsDose}];
        end
        
        %Set only the parameters describing some kind of reference dose,
        %where doseParams is an array of numeric values
        function obj = setDoseParameters(obj,doseParams)
            c = mat2cell(doseParams,1,numel(doseParams));
            [obj.parameters{obj.parameterIsDose}] = deal(c{:});
        end
    end
end

