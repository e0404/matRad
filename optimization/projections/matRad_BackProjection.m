classdef matRad_BackProjection
   
    properties (Access = private)
        wCache
        d
    end
    
    properties 
        dij
    end

    
    methods
        function obj = matRad_BackProjection()
            obj.wCache = [];
            obj.d = [];            
        end
        
        function obj = compute(obj,dij,w)
            if ~isequal(obj.wCache,w)
                obj.d = obj.computeResult(dij,w);
                obj.wCache = w;
            end
        end
        
        function d = GetResult(obj)
            d = obj.d;
        end
    end
    
    methods (Abstract)
        d = computeResult(obj,dij,w)
    end
end

