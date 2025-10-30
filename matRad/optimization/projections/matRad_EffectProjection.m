classdef matRad_EffectProjection < matRad_BackProjection
% matRad_EffectProjection class for effect-based optimization
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
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
            if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')
                if isempty(dij.mAlphaDose{scen}) || isempty(dij.mSqrtBetaDose{scen})
                    effect = [];
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispWarning('Empty mAlphaDose scenario in optimization detected! This should not happen...\n');
                else
                    effect = dij.mAlphaDose{scen}*w + (dij.mSqrtBetaDose{scen}*w).^2;
                end 
            else
                [ctScen,~,~] = ind2sub(size(dij.physicalDose),scen); %TODO: Workaround for now
                if isempty(dij.ax{ctScen}) || isempty(dij.bx{ctScen})
                    effect = [];
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispWarning('Empty dij.ax scenario in optimization detected! This should not happen...\n');
                elseif isfield(dij,'RBE') && (isscalar(dij.RBE) || numel(dij.RBE) == size(dij.physicalDose{scen},1)) && all(isfinite(dij.RBE))
                    effect = dij.ax{ctScen} .* (dij.physicalDose{scen} * w .* dij.RBE) + dij.bx{ctScen} .* (dij.physicalDose{scen} * w .* dij.RBE).^2;
                else
                    effect = dij.ax{ctScen} .* (dij.physicalDose{scen} * w) + dij.bx{ctScen} .* (dij.physicalDose{scen}*w).^2;
                end
            end
        end
        
        function wGrad = projectSingleScenarioGradient(~,dij,doseGrad,scen,w)
            if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')
                if isempty(dij.mAlphaDose{scen}) || isempty(dij.mSqrtBetaDose{scen})
                    wGrad = [];
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispWarning('Empty mAlphaDose scenario in optimization detected! This should not happen...\n');
                else
                    alphaTerm = (doseGrad{scen}' * dij.mAlphaDose{scen})';
                    quadTerm = dij.mSqrtBetaDose{scen} * w;
                    betaTerm = (2*(doseGrad{scen}.*quadTerm)' * dij.mSqrtBetaDose{scen})';
                    wGrad = alphaTerm + betaTerm;
                end
            else
                [ctScen,~,~] = ind2sub(size(dij.physicalDose),scen); %TODO: Workaround for now
                if isempty(dij.ax{ctScen}) || isempty(dij.bx{ctScen})
                    wGrad = [];
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispWarning('Empty dij.ax/dij.bx scenario in optimization detected! This should not happen...\n');
                else
                    physDose = dij.physicalDose{scen} * w;
                    if isfield(dij,'RBE') && (isscalar(dij.RBE) || numel(dij.RBE) == size(dij.physicalDose{scen},1)) && all(isfinite(dij.RBE))
                        alpha = dij.ax{ctScen} .* dij.RBE;
                        beta = dij.bx{ctScen} .* dij.RBE.^2;
                    else
                        alpha = dij.ax{ctScen};
                        beta = dij.bx{ctScen};
                    end
                    alphaTerm = ((doseGrad{scen} .* alpha)' * dij.physicalDose{scen})';
                    betaTerm = (2 * (doseGrad{scen} .* physDose .* beta)' * dij.physicalDose{scen})';
                    
                    wGrad = alphaTerm + betaTerm;
                end
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
                alphaTerm = (dExpGrad{scen}' * dij.mAlphaDoseExp{scen})';
                quadTerm = dij.mSqrtBetaDoseExp{scen} * w;
                betaTerm = (2*(dExpGrad{scen}.*quadTerm)' * dij.mSqrtBetaDoseExp{scen})';
                wGrad = alphaTerm + betaTerm;
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
