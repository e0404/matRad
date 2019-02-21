classdef matRad_OptimizationProblemDAO < matRad_OptimizationProblem
    %handle class to keep state easily
    
    %properties (Access = private)
    %    currentDose
    %    currentWeights
    %end
    properties
        apertureInfo
    end    
    
    methods (Static)
        %In External Files
        updatedInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVect);
        
        [apertureInfoVec, mappingMx, limMx] = matRad_daoApertureInfo2Vec(apertureInfo);
    end
    
    methods
        function obj = matRad_OptimizationProblemDAO(backProjection,apertureInfo)
            obj = obj@matRad_OptimizationProblem(backProjection);
            obj.apertureInfo = apertureInfo;
        end       
        
        function lb = lowerBounds(obj,w)
            lb = obj.apertureInfo.limMx(:,1);  % Lower bound on the variables. 
        end
        
        function ub = upperBounds(obj,w)            
            ub = obj.apertureInfo.limMx(:,2);
        end
    end
end

