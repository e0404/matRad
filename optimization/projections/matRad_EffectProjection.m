classdef matRad_EffectProjection < matRad_BackProjection
% matRad_EffectProjection class for effect-based optimization
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
        function obj = matRad_EffectProjection()
        end
        
        function e = computeResult(obj,dij,w)
            e = cell(numel(dij.mAlphaDose));
            % calculate effect
            for i = 1:numel(e)                
                e{i} = obj.computeSingleScenarioEffect(dij,w,i);                  
            end
        end
        
        function effect = computeSingleScenarioEffect(obj,dij,w,i)
            linTerm = dij.mAlphaDose{i} * w;
            quadTerm = dij.mSqrtBetaDose{i} * w;
            effect = linTerm + quadTerm.^2;
        end
    end
end
