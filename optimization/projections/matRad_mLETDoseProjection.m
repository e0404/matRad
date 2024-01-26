classdef matRad_mLETDoseProjection < matRad_BackProjection
% matRad_mLETProjection class to compute physical dose during optimization
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
        function obj = matRad_mLETDoseProjection()
            
        end
    end
    
    methods 
        function mLETd = computeSingleScenario(~,dij,scen,w)
            if ~isempty(dij.mLETDose{scen})
                mLETd = dij.mLETDose{scen}*w;
            else
                mLETd = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end 
        end
        
        function [dLDExp,dLDOmegaV] = computeSingleScenarioProb(~,dij,scen,w)
            if ~isempty(dij.mLETDoseExp{scen})
                dLDExp = dij.mLETDoseExp{scen}*w;

                for i = 1:size(dij.mLETDoseOmega,2)
                   dLDOmegaV{scen,i} = dij.mLETDoseOmega{scen,i} * w;
                end 
            else
                dLDExp = [];
                dLDOmegaV = [];
            end             
        end
     
        function wGrad = projectSingleScenarioGradient(~,dij,mLETDoseGrad,scen,~)
            if ~isempty(dij.mLETDose{scen})
               wGrad = (mLETDoseGrad{scen}' * dij.mLETDose{scen})';
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
        
        function wGrad = projectSingleScenarioGradientProb(~,dij,dLDExpGrad,dLDOmegaVgrad,scen,~)
            if ~isempty(dij.mLETDoseExp{scen})
                wGrad = (dLDExpGrad{scen}' * dij.mLETDoseExp{scen})';
                wGrad = wGrad + 2 * dLDOmegaVgrad;
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
    end
end


