classdef matRad_VariableRBEProjection < matRad_EffectProjection
% matRad_VariableRBEProjection class for RBE-weighted dose optimization
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
        function obj = matRad_VariableRBEProjection()
        end
        
        function RBExD = computeSingleScenario(obj,dij,scen,w)
            effect = computeSingleScenario@matRad_EffectProjection(obj,dij,scen,w); %First compute effect
            RBExD = zeros(dij.doseGrid.numOfVoxels,1);
            RBExD(dij.ixDose) = sqrt((effect(dij.ixDose)./dij.bx(dij.ixDose))+(dij.gamma(dij.ixDose).^2)) - dij.gamma(dij.ixDose);
        end
        
        function wGrad = projectSingleScenarioGradient(obj,dij,doseGrad,scen,w)
            if isempty(dij.mAlphaDose{scen}) || isempty(dij.mSqrtBetaDose{scen})
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            else
                %While the dose cache should be up to date here, we ask for
                %a computation (will skip if weights are equal to cache)
                obj = obj.compute(dij,w);
                
                %Scaling vor variable RBExD
                scaledEffect = obj.d{scen} + dij.gamma;
                doseGradTmp = zeros(dij.doseGrid.numOfVoxels,1);
                doseGradTmp(dij.ixDose) = doseGrad{scen}(dij.ixDose) ./ (2*dij.bx(dij.ixDose).*scaledEffect(dij.ixDose));
                
                %Now modify the effect computation
                vBias = (doseGradTmp' * dij.mAlphaDose{scen})';
                quadTerm = dij.mSqrtBetaDose{scen} * w;
                mPsi = (2*(doseGrad{scen}.*quadTerm)' * dij.mSqrtBetaDose{scen})';
                wGrad = vBias + mPsi;
            end
        end
        
        function [RBExDexp,dOmegaV] = computeSingleScenarioProb(~,dij,scen,w)
            if isempty(dij.mAlphaDoseExp{scen}) || isempty(dij.mSqrtBetaDoseExp{scen})
                RBExDexp = [];
                dOmegaV = [];
                %matRad_cfg = MatRad_Config.instance();
                %matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Probabilistic Backprojection uses inaccurate approximation for Variable RBE computation...\n');
                eExpLinTerm = dij.mAlphaDose{scen}*w;
                eExpSqTerm  = dij.mSqrtBetaDose{scen}*w;
                eExp = eExpLinTerm + eExpSqTerm.^2;
                
                RBExDexp = zeros(dij.doseGrid.numOfVoxels,1);
                RBExDexp(dij.ixDose) = sqrt((eExp(dij.ixDose)./dij.bx(dij.ixDose))+(dij.gamma(dij.ixDose).^2)) - dij.gamma(dij.ixDose);
                
                for i = 1:size(dij.physicalDoseOmega,2)
                   dOmegaV{scen,i} = dij.mAlphaDoseOmega{scen,i} * w;
                end 
            end      
        end
        
        function wGrad = projectSingleScenarioGradientProb(obj,dij,dExpGrad,dOmegaVgrad,scen,~)
            if isempty(dij.mAlphaDoseExp{scen}) || isempty(dij.mSqrtBetaDoseExp{scen})
                wGrad = [];
            else
                %While the dose cache should be up to date here, we ask for
                %a computation (will skip if weights are equal to cache)
                obj = obj.compute(dij,w);
                
                %Scaling vor variable RBExD
                scaledEffect = obj.dExp{scen} + dij.gamma;
                doseGradTmp = zeros(dij.doseGrid.numOfVoxels,1);
                doseGradTmp(dij.ixDose) = dExpGrad{scen}(dij.ixDose) ./ (2*dij.bx(dij.ixDose).*scaledEffect(dij.ixDose));
                
                %Now modify the effect computation
                vBias = (doseGradTmp' * dij.mAlphaDoseExp{scen})';
                quadTerm = dij.mSqrtBetaDose{scen} * w;
                mPsi = (2*(doseGradTmp.*quadTerm)' * dij.mSqrtBetaDoseExp{scen})';
                wGrad = vBias + mPsi + 2 * dOmegaVgrad;
            end
        end
    end
    
    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc,~,~)
            %Do nothing here to overwrite the behavior of the Effect
            %projection, since we have unit GyRBE here
        end
    end
end

