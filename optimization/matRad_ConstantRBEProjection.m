classdef matRad_ConstantRBEProjection < matRad_BackProjection
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
        
    methods
        function obj = matRad_ConstantRBEProjection()
        end
        
        function d = computeResult(obj,w,dij)
            d = cell(numel(dij.physicalDose));
            for i = numel(dij.physicalDose)                
                d{i} =  dij.physicalDose{i} * (w * dij.RBE);                
            end
        end
    end
end

