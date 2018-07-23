classdef (Abstract) matRad_Optimizer
    %OPTIMIZER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract)
        options %options struct
    end
    
    methods
        function obj = matRad_Optimizer
            %OPTIMIZER Construct an instance of this class
            %   Detailed explanation goes here
            obj = createDefaultOptimizerOptions(obj);
        end
        
        %function result = optimize(
    end
    
    methods (Abstract)
        obj  = createDefaultOptimizerOptions(obj)
        
        
    end
end

