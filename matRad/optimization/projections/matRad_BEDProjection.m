classdef matRad_BEDProjection < matRad_EffectProjection
% matRad_BEDProjection class for BED-based optimization
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
        function obj = matRad_BEDProjection()
        end
    end

    methods
        function BED = computeSingleScenario(obj,dij,scen,w)
                effect = computeSingleScenario@matRad_EffectProjection(obj,dij,scen,w);
                BED = zeros(dij.doseGrid.numOfVoxels,1);
                ix = ~(dij.ax{scen}==0);
                BED(ix) = effect(ix)./ dij.ax{scen}(ix);
                % photon equivalent BED = n * effect / alphax
            end 

        function wGrad = projectSingleScenarioGradient(obj,dij,doseGrad,scen,w)
                ix = ~(dij.ax{scen}==0);
                doseGradtmp{scen} = zeros(size(doseGrad{scen}));
                doseGradtmp{scen}(ix) = doseGrad{scen}(ix)./dij.ax{scen}(ix);
                wGradEffect = projectSingleScenarioGradient@matRad_EffectProjection(obj,dij,doseGradtmp,scen,w); 
                wGrad = wGradEffect;
            end
        end
end
