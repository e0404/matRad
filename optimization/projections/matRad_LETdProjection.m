classdef matRad_LETdProjection < matRad_BackProjection
% matRad_LETdProjection class to compute physical dose during optimization
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
        function obj = matRad_LETdProjection()
            
        end
    end
    
    methods 
        function LETd = computeSingleScenario(~,dij,scen,wCombined)

             if ~isempty(dij.mLETDose{scen})
                % separation of combined quantities
                wFinal = size(dij.physicalDose{scen},2); % how many pencil beams I have;
                w = wCombined(1:wFinal);
                I = wCombined(wFinal+1:end);
                d = dij.physicalDose{scen}*w;
                
                % Computation of the LETd
                LETD = dij.mLETDose{scen}.^2 * (w.*I); % computes the nominator quickly (should be correct)
                LETD(d > 0,1) = LETD(d > 0)./d(d > 0); %dose averging -> avoid div by zero
              
                LETd = [d;LETD];
                
            else
                LETd = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end

        end

        function [lExp,lOmegaV] = computeSingleScenarioProb(~,dij,scen,w)
            if ~isempty(dij.LETdExp{scen})
                lExp = dij.LETdExp{scen}*w;

                for i = 1:size(dij.LETdOmega,2)
                   lOmegaV{scen,i} = dij.LETdOmega{scen,i} * w;
                end 
            else
                lExp = [];
                lOmegaV = [];
            end             
        end

        function wGrad = projectSingleScenarioGradient(this,dij,doseGradCombined,scen,wCombined)
             if ~isempty(dij.mLETDose{scen})
               
                LETd = this.computeSingleScenario(dij,scen,wCombined);
                
                wFinal = size(dij.physicalDose{scen},2); % how many pencil beams I have;
                w = wCombined(1:wFinal);
                I = wCombined(1+wFinal:end);
                
                %doseGrad = doseGradCombined{scen}(1:dij.doseGrid.numOfVoxels,:);
                LETdGrad = doseGradCombined{scen}(dij.doseGrid.numOfVoxels+1:end,:);
                
                dose = LETd(1:dij.doseGrid.numOfVoxels);
                doseRate = LETd(dij.doseGrid.numOfVoxels+1:end);
                ix = dose > 0;
                tmp = zeros(size(dose));
                tmp(ix) = LETdGrad(ix) ./ dose(ix);
                tmpForW = zeros(size(tmp));
                tmpForW(ix) = tmp(ix) ./ dose(ix);

                %wGrad= (doseGrad' * dij.physicalDose{scen})' + LETdGradPart.*w;
                wGradPart1 = ((tmpForW.*dose)'  * dij.mLETDose{scen}.^2)' .* I; %u'v/v^2
                wGradPart2 = ((tmpForW.*doseRate)' * dij.mLETDose{scen})'; %uv'/v^2 

                wGradient = wGradPart1 - wGradPart2;

                LETdGradPart = (tmp' * dij.mLETDose{scen}.^2)';                                
                IGrad = LETdGradPart;               
                
                wGrad = [wGradient; IGrad];
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
            
        end

        function wGrad = projectSingleScenarioGradientProb(~,dij,lExpGrad,lOmegaVgrad,scen,~)
            if ~isempty(dij.LETdExp{scen})
                wGrad = (lExpGrad{scen}' * dij.LETdExp{scen})';
                wGrad = wGrad + 2 * lOmegaVgrad;
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
    end
end