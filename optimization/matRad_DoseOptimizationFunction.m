classdef (Abstract) matRad_DoseOptimizationFunction
    %MATRAD_DOSEOPTIMIZATIONFUNCTION This is the superclass for all
    %objectives and constraints to enable easy one-line identification
        
    properties (Abstract, Constant)        
        name                %Display name of the Objective. Needs to be implemented in sub-classes.
        parameterNames      %Cell array of Display names of the parameters. Needs to be implemented in sub-classes.
        parameterTypes      %Cell array of parameter types. Valid types are 'dose', 'numeric', or a cell list of string options. Needs to be implemented in sub-classes.
    end
    
    properties (Abstract, Access = public)
        parameters  
    end
    
    methods
        function obj = matRad_DoseOptimizationFunction(dataStruct)
            if isstruct(dataStruct)
                for fn = fieldnames(dataStruct)'    %enumerat fields
                    try
                        object.(fn{1}) = dataStruct.(fn{1});   %and copy
                    catch
                        warning('Could not copy field %s', fn{1});
                    end
                end
            end
        end
    end
    
    
    %Helper methods
    methods (Access = public)
        function doseParams = getDoseParameters(obj)
            %Get only the dose related parameters.
            ix = cellfun(@(c) isequal('dose',c),obj.parameterTypes);
            doseParams = [obj.parameters{ix}];
        end
        
        function obj = setDoseParameters(obj,doseParams)
            %Set only the dose related parameters.
            ix = cellfun(@(c) isequal('dose',c),obj.parameterTypes);
            obj.parameters(ix) = num2cell(doseParams);

        end
    end
end

