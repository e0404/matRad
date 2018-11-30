classdef (Abstract) matRad_Optimizer < handle
    %OPTIMIZER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract)
        options %options struct
        wResult
        resultInfo
    end
    
    methods
        function obj = matRad_Optimizer
            %OPTIMIZER Construct an instance of this class
            %   Detailed explanation goes here
            obj = createDefaultOptimizerOptions(obj);
        end
        
        
    end
    
    methods (Abstract)
        obj  = createDefaultOptimizerOptions(obj)
        
        obj = optimize(obj,w0,optiProb,dij,cst)
        
        [msg,statusflag] = GetStatus(obj)
    end
end

