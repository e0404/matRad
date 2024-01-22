classdef matRad_LETtProjection < matRad_BackProjection
% matRad_LETtProjection class to compute track averaged LET during optimization
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
        function obj = matRad_LETtProjection()
            
        end
    end
    
    methods 
        function LETt = computeSingleScenario(~,dij,scen,w)

            if ~isempty(dij.mLETDose{scen})

                LETw = dij.LETvD * w;
                log = logical(dij.LETvD) * w;
                LETt = LETw ./ log;

            else
                LETt = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end

        end

        function [ltExp,ltOmegaV] = computeSingleScenarioProb(~,dij,scen,w)
            if ~isempty(dij.LETtExp{scen})
                ltExp = dij.LETtExp{scen}*w;

                for i = 1:size(dij.LETtOmega,2)
                   ltOmegaV{scen,i} = dij.LETtOmega{scen,i} * w;
                end 
            else
                ltExp = [];
                ltOmegaV = [];
            end             
        end

        function wGrad = projectSingleScenarioGradient(~,dij,LETtGrad,scen,w)
             if ~isempty(dij.mLETDose{scen})

                log = logical(dij.LETvD) * w;
                log2 = log.^2;

                vox_tmp = zeros(dij.doseGrid.numOfVoxels,1);
                vox_tmp =  dij.LETvD * w;
                firstDerivativeterm = (vox_tmp .* LETtGrad{1}) ./ log2;
                
                vox_tmp = zeros(dij.doseGrid.numOfVoxels,1);
                vox_tmp = dij.LETvD .* log;
                secondDerivativeterm =  ((vox_tmp .* LETtGrad{1}) * w) ./ log2; % there shouldn't be an extra w but otherwise the sizes don't match
               
                wGrad = (firstDerivativeterm - secondDerivativeterm);
              
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
            
        end

        function wGrad = projectSingleScenarioGradientProb(~,dij,ltExpGrad,ltOmegaVgrad,scen,~)
            if ~isempty(dij.LETtExp{scen})
                wGrad = (ltExpGrad{scen}' * dij.LETtExp{scen})';
                wGrad = wGrad + 2 * ltOmegaVgrad;
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
    end
end