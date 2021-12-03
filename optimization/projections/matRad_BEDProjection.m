classdef matRad_BEDProjection < matRad_EffectProjection
% matRad_BEDProjection class for BED-based optimization
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    methods
        function obj = matRad_BEDProjection()
        end
    end
    
    methods
        function BED = computeSingleScenario(obj,dij,scen,w)
                effect = computeSingleScenario@matRad_EffectProjection(obj,dij,scen,w);
                BED = zeros(dij.doseGrid.numOfVoxels,1);
                ix = ~(dij.ax==0);
                BED(ix) = obj.numOfFractions.*effect(ix)./ dij.ax(ix);
                % photon equivalent BED = n * effect / alphax
            end 
        
        function wGrad = projectSingleScenarioGradient(obj,dij,doseGrad,scen,w)
                doseGrad{scen} = doseGrad{scen}./dij.ax;
                wGradEffect = projectSingleScenarioGradient@matRad_EffectProjection(obj,dij,doseGrad,scen,w); 
                wGrad = obj.numOfFractions.*wGradEffect;
            end
        end
end

