classdef matRad_VariableRBEProjection < matRad_EffectProjection
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
        
    methods
        function obj = matRad_VariableRBEProjection()
        end
        
        function RBExD = computeResult(obj,dij,w)
            %effect = computeResult@matRad_EffectProjection(obj,dij,w)
            RBExD = cell(size(dij.mAlphaDose));
            % calculate effect
            for i = 1:numel(RBExD)                 
                e = obj.computeSingleScenarioEffect(dij,w,i);
                RBExD{i} = zeros(dij.numOfVoxels,1);
                RBExD{i}(dij.ixDose) = sqrt((e(dij.ixDose)./dij.bx(dij.ixDose))+(dij.gamma(dij.ixDose).^2)) - dij.gamma(dij.ixDose);
            end
        end
    end
end

