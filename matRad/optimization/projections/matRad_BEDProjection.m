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
            %Get corresponding ct scenario
            [ctScen,~] = ind2sub(size(dij.physicalDose),scen);

            effect = computeSingleScenario@matRad_EffectProjection(obj,dij,scen,w);
            
            BED = zeros(dij.doseGrid.numOfVoxels,1);
            BED(dij.ixDose{ctScen}) = effect(dij.ixDose{ctScen})./ dij.ax{ctScen}(dij.ixDose{ctScen});
            % photon equivalent BED = n * effect / alphax
        end

        function wGrad = projectSingleScenarioGradient(obj,dij,doseGrad,scen,w)
            [ctScen,~] = ind2sub(size(dij.physicalDose),scen);

            doseGradtmp{scen} = zeros(size(doseGrad{scen}));
            doseGradtmp{scen}(dij.ixDose{ctScen}) = doseGrad{scen}(dij.ixDose{ctScen})./dij.ax{scen}(dij.ixDose{ctScen});
            wGradEffect = projectSingleScenarioGradient@matRad_EffectProjection(obj,dij,doseGradtmp,scen,w);
            wGrad = wGradEffect;
        end
    end
    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc,alphaX,betaX)
            doses = optiFunc.getDoseParameters();    
            BED = doses*(1 + doses/(alphaX/betaX));        
            optiFunc = optiFunc.setDoseParameters(BED);
        end
    end
end
