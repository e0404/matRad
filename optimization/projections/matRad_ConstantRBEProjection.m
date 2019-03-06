classdef matRad_ConstantRBEProjection < matRad_BackProjection
% matRad_BackProjection for optimization based on RBExDose
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
        function obj = matRad_ConstantRBEProjection()
        end
        
        function d = computeResult(obj,dij,w)
            d = cell(numel(dij.physicalDose));
            for i = numel(dij.physicalDose)                
                d{i} =  dij.physicalDose{i} * (w * dij.RBE);                
            end
        end
    end
end

