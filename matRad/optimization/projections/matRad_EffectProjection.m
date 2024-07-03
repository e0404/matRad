classdef matRad_EffectProjection < matRad_BackProjection
% matRad_EffectProjection class for effect-based optimization
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
        function obj = matRad_EffectProjection()
        end
    end
    
    methods
        function effect = computeSingleScenario(~,dij,scen,w)
            if isempty(dij.mAlphaDose{scen}) || isempty(dij.mSqrtBetaDose{scen})
                effect = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            else
                effect = dij.mAlphaDose{scen}*w + (dij.mSqrtBetaDose{scen}*w).^2;
            end 
        end
        
        function wGrad = projectSingleScenarioGradient(~,dij,doseGrad,scen,w)
            if isempty(dij.mAlphaDose{scen}) || isempty(dij.mSqrtBetaDose{scen})
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            else
                vBias = (doseGrad{scen}' * dij.mAlphaDose{scen})';
                quadTerm = dij.mSqrtBetaDose{scen} * w;
                mPsi = (2*(doseGrad{scen}.*quadTerm)' * dij.mSqrtBetaDose{scen})';
                wGrad = vBias + mPsi;
            end
        end
        
        function [eExp,dOmegaV] = computeSingleScenarioProb(~,dij,scen,w)
            if isempty(dij.mAlphaDoseExp{scen}) || isempty(dij.mSqrtBetaDoseExp{scen})
                eExp = [];
                dOmegaV = [];
                %matRad_cfg = MatRad_Config.instance();
                %matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Probabilistic Backprojection uses inaccurate approximation for effect computation...\n');
                %eExpLinTerm = dij.mAlphaDose{scen}*w;
                %eExpSqTerm  = dij.mSqrtBetaDose{scen}*w;
                eExp = dij.mAlphaDose{scen}*w + (dij.mSqrtBetaDose{scen}*w).^2;
                
                for i = 1:size(dij.physicalDoseOmega,2)
                   dOmegaV{scen,i} = dij.mAlphaDoseOmega{scen,i} * w;
                end 
            end      
        end
        
        function wGrad = projectSingleScenarioGradientProb(~,dij,dExpGrad,dOmegaVgrad,scen,~)
            if isempty(dij.mAlphaDoseExp{scen}) || isempty(dij.mSqrtBetaDoseExp{scen})
                wGrad = [];
            else
                vBias = (dExpGrad{scen}' * dij.mAlphaDoseExp{scen})';
                quadTerm = dij.mSqrtBetaDoseExp{scen} * w;
                mPsi = (2*(dExpGrad{scen}.*quadTerm)' * dij.mSqrtBetaDoseExp{scen})';
                wGrad = vBias + mPsi;
                wGrad = wGrad + 2 * dOmegaVgrad;
            end            
        end        
    end
    
    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc,alphaX,betaX)
            doses = optiFunc.getDoseParameters();    
            effect = alphaX*doses + betaX*doses.^2;        
            optiFunc = optiFunc.setDoseParameters(effect);
        end
    end
end
