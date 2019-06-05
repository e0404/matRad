classdef matRad_BackProjection
% matRad_BackProjection superclass for all backprojection algorithms used
% within matRad optimzation processes
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
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
    
    %These should be abstract methods, however Octave can't parse them. As soon 
    %as Octave is able to do this, they should be made abstract again 
    methods %(Abstract)
        function d = computeResult(obj,dij,w)
          error('Function needs to be implemented!');
        end
    end
end

