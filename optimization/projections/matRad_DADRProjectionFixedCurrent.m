classdef matRad_DADRProjectionFixedCurrent < matRad_BackProjection
    
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
        function obj = matRad_DADRProjectionFixedCurrent()
            
        end
    end
    
    methods
        function DADR = computeSingleScenario(this,dij,scen,w) 
            if ~isempty(dij.physicalDose{scen})
                I = dij.fixedCurrent;                
                nA2particles = 1e-9 ./ 1.602176634e-19 ./ 1e6;
                I = nA2particles * I;
                
                d = dij.physicalDose{scen}*w;
                %
                % Computation of the DADR
                DADR = dij.physicalDose{scen}.^2 * (I*w); %computes the nominator quickly (should be correct)
                DADR(d > 0,1) = DADR(d > 0)./d(d > 0); %dose averging -> avoid div by zero                                             
            else
                dCombined = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
        
        function wGrad = projectSingleScenarioGradient(this,dij,dadrGrad,scen,w)
            if ~isempty(dij.physicalDose{scen})
                %scalingFac = 10^4;
                
                %dadr = this.computeSingleScenario(dij,scen,w);
                dose = dij.physicalDose{scen}*w;
                %DADR = dij.physicalDose{scen}.^2 * (w.*I); %computes the nominator quickly (should be correct)
                %DADR(dose > 0,1) = DADR(dose > 0)./dose(dose > 0); %dose averging -> avoid div by zero  
                
                I = dij.fixedCurrent;
                nA2particles = 1e-9 ./ 1.602176634e-19 ./ 1e6;
                I = nA2particles * I;
                
                %Gradient computation
                ix = dose > 0;
                tmp = zeros(size(dose));
                tmp(ix) = dadrGrad{1}(ix) ./ dose(ix).^2;

                wGrad = tmp' * (dij.physicalDose{scen}.^2 .* dose - (dij.physicalDose{scen}.^2*w).*dij.physicalDose{scen});

                wGrad = I*wGrad';              

                %{
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
                %}
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
    end
    
    
end

