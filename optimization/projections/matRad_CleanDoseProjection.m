classdef matRad_CleanDoseProjection < matRad_BackProjection
% matRad_CleanDoseProjection class to compute clean dose during optimization
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
        function obj = matRad_CleanDoseProjection()
            
        end
    end
    
    methods 
        function d = computeSingleScenario(~,dij,scen,w)
            if ~isempty(dij.cleanDose{scen})
                d = dij.cleanDose{scen}*w;
            else
                d = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end 
        end
        
        function [cExp,cOmegaV] = computeSingleScenarioProb(~,dij,scen,w)
            if ~isempty(dij.cleanDoseExp{scen})
                cExp = dij.cleanDoseExp{scen}*w;

                for i = 1:size(dij.cleanDoseOmega,2)
                   cOmegaV{scen,i} = dij.cleanDoseOmega{scen,i} * w;
                end 
            else
                cExp = [];
                cOmegaV = [];
            end             
        end
     
        function wGrad = projectSingleScenarioGradient(~,dij,cleanDoseGrad,scen,~)
            if ~isempty(dij.cleanDose{scen})
               wGrad = (cleanDoseGrad{scen}' * dij.cleanDose{scen})';
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
        
        function wGrad = projectSingleScenarioGradientProb(~,dij,cExpGrad,cOmegaVgrad,scen,~)
            if ~isempty(dij.cleanDoseExp{scen})
                wGrad = (cExpGrad{scen}' * dij.cleanDoseExp{scen})';
                wGrad = wGrad + 2 * cOmegaVgrad;
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
    end
end


