classdef (Abstract) matRad_DoseConstraint
    %MATRAD_DOSEOBJECTIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract, Constant)        
        name
        parameterNames
        %parameterIsDose
        parameterTypes
    end
    
    properties (Abstract, Access = public)
        parameters
    end
        
    methods (Abstract)
        cDose        = computeDoseConstraintFunction(obj,dose)
        cDoseJacob   = computeDoseConstraintJacobian(obj,dose)
        cu           = upperBounds(obj,n)
        cl           = lowerBounds(obj,n)                
    end
    
        %Helper methods
    methods (Access = public)
        function jStruct = getDoseConstraintJacobianStructure(obj,n)
            %The default structure for scalar constraint functions is a
            %gradient vector
            jStruct = ones(n,1);
        end
  
        
        %Get only the parameters describing some kind of reference dose as
        %numeric array
        function doseParams = getDoseParameters(obj)
            ix = cellfun(@(c) isequal('dose',c),obj.parameterTypes);
            doseParams = [obj.parameters{ix}];
        end
        
        %Set only the parameters describing some kind of reference dose,
        %where doseParams is an array of numeric values
        function obj = setDoseParameters(obj,doseParams)
            %c = mat2cell(doseParams,1,numel(doseParams));
            %[obj.parameters{obj.parameterIsDose}] = deal(c{:});
            ix = cellfun(@(c) isequal('dose',c),obj.parameterTypes);
            obj.parameters(ix) = num2cell(doseParams);

        end
    end
end

