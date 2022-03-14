classdef matRad_XBDDADRProjection < matRad_BackProjection
% matRad_BackProjection for optimization based on RBExDose
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
        k = 0.5;
        DADR_t = 40;
        a_num = 4; %a = a_num / DADR_t
    end

    methods
        function obj = matRad_XBDDADRProjection()
        end   
    end
    
    methods 
        function XBD = computeSingleScenario(obj,dij,scen,w)
            if ~isempty(dij.physicalDose{scen})
                I = dij.fixedCurrent;                
                nA2particles = 1e-9 ./ 1.602176634e-19 ./ 1e6;
                I = nA2particles * I;
                
                d = dij.physicalDose{scen}*w;
                %
                % Computation of the DADR
                DADR = dij.physicalDose{scen}.^2 * (I*w); %computes the nominator quickly (should be correct)
                DADR(d > 0) = DADR(d > 0)./d(d > 0); %dose averging -> avoid div by zero 
                
                
                a = obj.a_num/obj.DADR_t;                
                XBD = d.*(1 - obj.k./(1+exp(-a*(DADR - obj.DADR_t))));
            else
                XBD = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
        
        function wGrad = projectSingleScenarioGradient(obj,dij,xbdGrad,scen,w)
            if ~isempty(dij.physicalDose{scen})
                %Current normalization
                I = dij.fixedCurrent;
                nA2particles = 1e-9 ./ 1.602176634e-19 ./ 1e6;
                I = nA2particles * I;

                %Required quantities
                dose = dij.physicalDose{scen}*w;
                dijSqW = dij.physicalDose{scen}.^2*w;
                DADR = dijSqW * I; 
                DADR(dose > 0) = DADR(dose > 0)./dose(dose > 0); %dose averging -> avoid div by zero 
                
                %Parameter a
                a = obj.a_num/obj.DADR_t;  
                
                %{
                %%THis tries to compute the derivative of XBD = d*(1-logistic(DADR))               
                %Product rule components
                %1) this should be the derivative of the dose part w.r.t. w multiplied with the logistic DADR function
                logTerm = 1 - obj.k./(1+exp(-a*(DADR - obj.DADR_t)));
                tmp = xbdGrad{scen}.*logTerm;
                duv = tmp'*(dij.physicalDose{scen}); 
                
                %2) This should compute the derivative of the logistic function 
                ix = dose > 0;
                tmp = tmp.*(dose.*(1 - logTerm));
                tmp(ix) = tmp(ix) * I ./ dose(ix).^2;    
                %uvd =tmp' * (dij.physicalDose{scen}.^2 .* dose - dij.physicalDose{scen} .* dijSqW);
                uvd = (tmp .* dose)' * dij.physicalDose{scen}.^2 - (tmp .* dijSqW)' * dij.physicalDose{scen};
                
                wGrad = (duv + uvd)';
                %wGrad = (duv + uvd)';                

                %u'v + uv'
                %wGrad = xbdGrad{scen}' * ( (logTerm .* dij.physicalDose{scen}) * logTermGradW )
                %}


                %%THis tries to compute the derivative of XBD = d - d*(logistic(DADR))  
                %Product rule components
                %1) this should be the derivative of the dose part w.r.t. w multiplied with the logistic DADR function
                logTerm = obj.k./(1+exp(-a*(DADR - obj.DADR_t)));
                tmp = xbdGrad{scen}.*logTerm;
                duv = tmp'*(dij.physicalDose{scen}); 
                
                %2) This should compute the derivative of the logistic function 
                ix = dose > 0;
                tmp = a*tmp.*(dose.*(1 - logTerm./obj.k));
                tmp(ix) = tmp(ix) * I ./ dose(ix).^2;    
                %uvd =tmp' * (dij.physicalDose{scen}.^2 .* dose - dij.physicalDose{scen} .* dijSqW);
                uvd = (tmp .* dose)' * dij.physicalDose{scen}.^2 - (tmp .* dijSqW)' * dij.physicalDose{scen};
                
                wGrad = (xbdGrad{scen}' * dij.physicalDose{scen} - duv - uvd)';

            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
    end
    

end

