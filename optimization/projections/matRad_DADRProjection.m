classdef matRad_DADRProjection < matRad_BackProjection
    
    % matRad_DoseProjection class to compute physical dose and DADR during optimization
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
    
    properties
    end


    methods
        function obj = matRad_DADRProjection()
            
        end
    end
    
    methods
        function dCombined = computeSingleScenario(this,dij,scen,wCombined) 
            if ~isempty(dij.physicalDose{scen})
                %separation of combined quantities
                wFinal = size(dij.physicalDose{scen},2); % how many pencil beams I have;
                w = wCombined(1:wFinal);
                I = wCombined(wFinal+1:end);
                d = dij.physicalDose{scen}*w;
                %
                % Computation of the DADR
                DADR = dij.physicalDose{scen}.^2 * (w.*I); %computes the nominator quickly (should be correct)
                DADR(d > 0,1) = DADR(d > 0)./d(d > 0); %dose averging -> avoid div by zero
              
                dCombined = [d;DADR];
                
            else
                dCombined = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
        
        function wTildaGrad = projectSingleScenarioGradient(this,dij,doseGradCombined,scen,wCombined)
            if ~isempty(dij.physicalDose{scen})
                %scalingFac = 10^4;
                
                dCombined = this.computeSingleScenario(dij,scen,wCombined);
                
                wFinal = size(dij.physicalDose{scen},2); % how many pencil beams I have;
                w = wCombined(1:wFinal);
                I = wCombined(1+wFinal:end);
                
                %doseGrad = doseGradCombined{scen}(1:dij.doseGrid.numOfVoxels,:);
                dadrGrad = doseGradCombined{scen}(dij.doseGrid.numOfVoxels+1:end,:);
                
                dose = dCombined(1:dij.doseGrid.numOfVoxels);
                doseRate = dCombined(dij.doseGrid.numOfVoxels+1:end);
                ix = dose > 0;
                tmp = zeros(size(dose));
                tmp(ix) = dadrGrad(ix) ./ dose(ix);
                tmpForW = zeros(size(tmp));
                tmpForW(ix) = tmp(ix) ./ dose(ix);

                %wGrad= (doseGrad' * dij.physicalDose{scen})' + dadrGradPart.*w;
                wGradPart1 = ((tmpForW.*dose)'  * dij.physicalDose{scen}.^2)' .* I; %u'v/v^2
                wGradPart2 = ((tmpForW.*doseRate)' * dij.physicalDose{scen})'; %uv'/v^2 

                wGrad = wGradPart1 - wGradPart2;

                dadrGradPart = (tmp' * dij.physicalDose{scen}.^2)';                                
                IGrad = dadrGradPart;               
                
                wTildaGrad = [wGrad; IGrad];
            else
                wTildaGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
    end
    
    
end

