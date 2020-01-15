classdef matRad_VariableRBEProjection < matRad_EffectProjection
% matRad_VariableRBEProjection class for RBE-weighted dose optimization
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
        function obj = matRad_VariableRBEProjection()
        end
        
        function RBExD = computeResult(obj,dij,w)
            %effect = computeResult@matRad_EffectProjection(obj,dij,w)
            RBExD = cell(size(dij.mAlphaDose));
            % calculate effect
            for i = 1:numel(RBExD)                 
                e = obj.computeSingleScenarioEffect(dij,w,i);
                RBExD{i} = zeros(dij.doseGrid.numOfVoxels,1);
                RBExD{i}(dij.ixDose) = sqrt((e(dij.ixDose)./dij.bx(dij.ixDose))+(dij.gamma(dij.ixDose).^2)) - dij.gamma(dij.ixDose);
            end
        end
    end
end

