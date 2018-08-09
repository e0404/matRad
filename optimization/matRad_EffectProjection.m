classdef matRad_EffectProjection < matRad_BackProjection
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
        
    methods
        function obj = matRad_EffectProjection()
        end
        
        function e = computeResult(obj,w,dij)
            e = cell(numel(dij.mAlphaDose));
            % calculate effect
            for i = 1:numel(e)                
                linTerm  = dij.mAlphaDose{i} * w;
                quadTerm = dij.mSqrtBetaDose{i} * w;
                e{i} = linTerm + quadTerm.^2;                   
            end
        end
    end
end
