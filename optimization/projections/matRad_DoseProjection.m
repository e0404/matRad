classdef matRad_DoseProjection < matRad_BackProjection
% matRad_DoseProjection class to compute physical dose during optimization
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
    
    methods
        function obj = matRad_DoseProjection()
            
        end
        
        function d = computeResult(obj,dij,w)
            d = cell(numel(dij.physicalDose));
            for i = numel(dij.physicalDose)
                d{i} = dij.physicalDose{i} * w;
            end
        end
    end
    

end

