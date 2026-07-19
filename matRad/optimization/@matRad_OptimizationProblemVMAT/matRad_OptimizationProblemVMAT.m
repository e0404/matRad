classdef matRad_OptimizationProblemVMAT < matRad_OptimizationProblemDAO
    % handle class to keep state easily

    methods (Static)
        % In External Files
        updatedInfo = matRad_daoVec2ApertureInfo(apertureInfo, apertureInfoVect)

        [apertureInfoVec, mappingMx, limMx] = matRad_daoApertureInfo2Vec(apertureInfo)

        % Aperture-info delivery post-processing (shared by the VMAT
        % sequencer and the DAO optimizer)
        apertureInfo = leafTouching(apertureInfo)

        apertureInfo = maxLeafSpeed(apertureInfo)

        apertureInfo = optDelivery(apertureInfo, fast)
    end

    methods

        function obj = matRad_OptimizationProblemVMAT(backProjection, apertureInfo)
            obj = obj@matRad_OptimizationProblemDAO(backProjection, apertureInfo);
        end

    end
end
