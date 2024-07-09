classdef matRad_LETxDoseProjection < matRad_BackProjection
% matRad_LETxProjection class to compute physical dose during optimization
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
        function obj = matRad_LETxDoseProjection()
            
        end
    end
    
    methods 
        function LxD = computeSingleScenario(~,dij,scen,w)
            if ~isempty(dij.mLETDose{scen})
                LxD = dij.mLETDose{scen}*w;
            else
                LxD = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end 
        end
        
        function [LxDExp,LxDOmegaV] = computeSingleScenarioProb(~,dij,scen,w)
            if ~isempty(dij.mLETDoseExp{scen})
                LxDExp = dij.mLETDoseExp{scen}*w;

                for i = 1:size(dij.mLETDoseOmega,2)
                   LxDOmegaV{scen,i} = dij.mLETDoseOmega{scen,i} * w;
                end 
            else
                LxDExp = [];
                LxDOmegaV = [];
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
        
        function wGrad = projectSingleScenarioGradientProb(~,dij,LxDExpGrad,LxDOmegaVgrad,scen,~)
            if ~isempty(dij.mLETDoseExp{scen})
                wGrad = (LxDExpGrad{scen}' * dij.mLETDoseExp{scen})';
                wGrad = wGrad + 2 * LxDOmegaVgrad;
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
    end
end


