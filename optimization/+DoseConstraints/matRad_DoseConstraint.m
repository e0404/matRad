classdef (Abstract) matRad_DoseConstraint
    %MATRAD_DOSEOBJECTIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract, Constant)        
        name
    end
    
    properties (Constant)
        type = 'Constraint';
    end
    
    properties (Abstract, Access = public)
        parameters        
    end
    
    methods (Abstract)
        cDose        = computeDoseConstraintFunction(obj,dose)
        cDoseJacob   = computeDoseConstraintJacobian(obj,dose)
        cu           = upperBounds(obj)
        cl           = lowerBounds(obj)
    end
    
end

