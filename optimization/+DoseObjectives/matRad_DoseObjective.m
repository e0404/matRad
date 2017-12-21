classdef (Abstract) matRad_DoseObjective
    %MATRAD_DOSEOBJECTIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract, Constant)        
        name
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
    
end

