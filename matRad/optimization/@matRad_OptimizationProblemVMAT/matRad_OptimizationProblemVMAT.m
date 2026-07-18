classdef matRad_OptimizationProblemVMAT < matRad_OptimizationProblemDAO
    % handle class to keep state easily

    methods (Static)
        % In External Files
        updatedInfo = matRad_daoVec2ApertureInfo(apertureInfo, apertureInfoVect)

        [apertureInfoVec, mappingMx, limMx] = matRad_daoApertureInfo2Vec(apertureInfo)
    end

    methods

        function obj = matRad_OptimizationProblemVMAT(backProjection, apertureInfo)
            obj = obj@matRad_OptimizationProblemDAO(backProjection, apertureInfo);
        end

    end
end
