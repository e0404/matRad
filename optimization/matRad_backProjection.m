classdef matRad_BackProjection
   
    properties (Access = private)
        w
        d
    end
    
    methods
        function obj = matRad_BackProjection()
            obj.w = [];
            obj.d = [];
            
        end
        
        function obj = compute(obj,w,dij)
            if ~isequal(obj.w,w)
                obj.d = obj.computeResult(w,dij);
                obj.w = w;
            end
        end
        
        function d = GetResult(obj)
            d = obj.d;
        end
    end
    
    methods (Abstract)
        d = computeResult(obj,w,dij)
    end
end

