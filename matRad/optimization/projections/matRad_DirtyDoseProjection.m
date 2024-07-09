classdef matRad_DirtyDoseProjection < matRad_BackProjection
% matRad_DirtyDoseProjection class to compute dirty dose during optimization
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
        function obj = matRad_DirtyDoseProjection()
            
        end
    end
    
    methods 
        function dD = computeSingleScenario(~,dij,scen,w)
            if ~isempty(dij.dirtyDose{scen})
                dD = dij.dirtyDose{scen}*w;
            else
                dD = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end 
        end
        
        function [dDExp,dDOmegaV] = computeSingleScenarioProb(~,dij,scen,w)
            if ~isempty(dij.dirtyDoseExp{scen})
                dDExp = dij.dirtyDoseExp{scen}*w;

                for i = 1:size(dij.dirtyDoseOmega,2)
                   dDOmegaV{scen,i} = dij.dirtyDoseOmega{scen,i} * w;
                end 
            else
                dDExp = [];
                dDOmegaV = [];
            end             
        end
     
        function wGrad = projectSingleScenarioGradient(~,dij,dirtyDoseGrad,scen,~)
            if ~isempty(dij.dirtyDose{scen})
               wGrad = (dirtyDoseGrad{scen}' * dij.dirtyDose{scen})';
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
        
        function wGrad = projectSingleScenarioGradientProb(~,dij,dDExpGrad,dDOmegaVgrad,scen,~)
            if ~isempty(dij.dirtyDoseExp{scen})
                wGrad = (dDExpGrad{scen}' * dij.dirtyDoseExp{scen})';
                wGrad = wGrad + 2 * dDOmegaVgrad;
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
    end
end


