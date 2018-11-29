classdef matRad_EffectProjection < matRad_BackProjection
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
        
    methods
        function obj = matRad_EffectProjection()
        end
        
        function e = computeResult(obj,dij,w)
            e = cell(numel(dij.mAlphaDose));
            % calculate effect
            for i = 1:numel(e)                
                e{i} = obj.computeSingleScenarioEffect(dij,w,i);                  
            end
        end
        
        function effect = computeSingleScenarioEffect(obj,dij,w,i)
            linTerm = dij.mAlphaDose{i} * w;
            quadTerm = dij.mSqrtBetaDose{i} * w;
            effect = linTerm + quadTerm.^2;
        end
    end
end
