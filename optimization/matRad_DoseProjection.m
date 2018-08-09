classdef matRad_DoseProjection < matRad_BackProjection
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = matRad_DoseProjection()
        end
        
        function d = computeResult(obj,w,dij)
            d = cell(numel(dij.physicalDose));
            for i = numel(dij.physicalDose)
                d{i} = dij.physicalDose{1} * w;
            end
        end
    end
    

end

