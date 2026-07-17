classdef matRad_DemonsImageRegistration < matRad_ImageRegistrationBase
    % matRad_DemonsImageRegistration Demons-based elastic image registration
    %   Computes deformation vector fields between the reference CT scenario
    %   and all other scenarios of a 4D CT using MATLAB's imregdemons
    %   (requires the Image Processing Toolbox) and can propagate contours
    %   from the reference scenario using push DVFs.
    %
    % call
    %   reg = matRad_DemonsImageRegistration()
    %   reg = matRad_DemonsImageRegistration(propStruct)
    %
    %   ct = reg.calcDVF(ct)
    %   [ct, cst] = reg.propContours(ct, cst)
    %
    % input
    %   propStruct: (optional) struct with registration settings, fields
    %               matching the class properties (e.g. from struct(reg))
    %
    % properties
    %   refScen:                    reference CT scenario index (default: 1)
    %   dvfType:                    'pull' or 'push' (default: 'pull')
    %   dvfUnits:                   'mm' or 'voxel' (default: 'mm')
    %   numIterations:              number of demons iterations (default: 100)
    %   pyramidLevels:              number of multi-resolution levels (default: 1)
    %   accumulatedFieldSmoothing:  smoothing of the accumulated field (default: 1)
    %
    % References
    %   https://www.mathworks.com/help/images/ref/imregdemons.html
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2026 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Constant)
        name = 'Demons Image Registration (imregdemons)'
    end

    properties
        numIterations = 100             % number of demons iterations per pyramid level
        pyramidLevels = 1               % number of multi-resolution pyramid levels
        accumulatedFieldSmoothing = 1   % gaussian smoothing of the accumulated displacement field
    end

    methods

        function this = matRad_DemonsImageRegistration(propStruct)
            if nargin < 1
                propStruct = [];
            end
            this@matRad_ImageRegistrationBase(propStruct);
        end

        function set.numIterations(this, value)
            matRad_cfg = MatRad_Config.instance();
            if ~(isnumeric(value) && isscalar(value) && isfinite(value) && round(value) == value && value >= 1)
                matRad_cfg.dispError('numIterations must be a positive integer scalar.');
            end
            this.numIterations = double(value);
        end

        function set.pyramidLevels(this, value)
            matRad_cfg = MatRad_Config.instance();
            if ~(isnumeric(value) && isscalar(value) && isfinite(value) && round(value) == value && value >= 1)
                matRad_cfg.dispError('pyramidLevels must be a positive integer scalar.');
            end
            this.pyramidLevels = double(value);
        end

        function set.accumulatedFieldSmoothing(this, value)
            matRad_cfg = MatRad_Config.instance();
            if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0)
                matRad_cfg.dispError('accumulatedFieldSmoothing must be a nonnegative scalar.');
            end
            this.accumulatedFieldSmoothing = double(value);
        end

        function ct = calcDVF(this, ct)
            matRad_cfg = MatRad_Config.instance();

            this.validateCt(ct);
            if strcmp(this.dvfUnits, 'mm')
                this.validateResolution(ct);
            end
            this.checkDependencies(false);

            ct.dvf = cell(1, ct.numOfCtScen);
            ct.dvfMetadata = struct('dvfType', this.dvfType, ...
                                    'dvfUnits', this.dvfUnits, ...
                                    'refScen', this.refScen);

            for scen = 1:ct.numOfCtScen
                if scen == this.refScen
                    ct.dvf{scen} = zeros([3 reshape(ct.cubeDim, 1, [])]);
                else
                    matRad_cfg.dispInfo('Registering CT scenario %d against reference scenario %d.\n', scen, this.refScen);
                    ct.dvf{scen} = this.calculateScenarioDvf(ct, scen);
                end
            end
        end

        function [ct, cst] = propContours(this, ct, cst)
            matRad_cfg = MatRad_Config.instance();

            this.validateCt(ct);

            if ~isfield(ct, 'dvf') || isempty(ct.dvf)
                if ~strcmp(this.dvfType, 'push')
                    matRad_cfg.dispError('Contour propagation requires push DVFs, but dvfType is ''%s''.', this.dvfType);
                end
                ct = this.calcDVF(ct);
            end

            metadata = this.getDvfMetadata(ct);
            if ~strcmp(metadata.dvfType, 'push')
                matRad_cfg.dispError('Contour propagation requires push DVFs, but the ct contains ''%s'' DVFs.', metadata.dvfType);
            end

            this.checkDependencies(true);
            this.validateDvfContainer(ct);
            refScenUsed = metadata.refScen;

            for structure = 1:size(cst, 1)
                if size(cst, 2) < 4 || isempty(cst{structure, 4}) || numel(cst{structure, 4}) < refScenUsed
                    matRad_cfg.dispError('Structure %d does not contain contours for reference CT scenario %d.', structure, refScenUsed);
                end

                refContour = cst{structure, 4}{1, refScenUsed};

                if isempty(refContour)
                    matRad_cfg.dispWarning('Structure %d has no reference contour. Leaving propagated contours empty.', structure);
                    for scen = 1:ct.numOfCtScen
                        if scen ~= refScenUsed
                            cst{structure, 4}{1, scen} = [];
                        end
                    end
                    continue
                end

                referenceMask = false(reshape(ct.cubeDim, 1, []));
                referenceMask(refContour) = true;
                matRad_cfg.dispInfo('Propagating contours of structure %d.\n', structure);

                for scen = 1:ct.numOfCtScen
                    if scen == refScenUsed
                        continue
                    end

                    % imwarp expects the displacement field in voxel units
                    % with size [cubeDim 3] and component order x,y,z
                    scenarioDvf = this.convertDvfToVoxelUnits(ct.dvf{scen}, ct, metadata.dvfUnits);
                    warpedMask = imwarp(single(referenceMask), permute(scenarioDvf, [2 3 4 1]));
                    cst{structure, 4}{1, scen} = find(warpedMask > 0.5);
                end
            end
        end

    end

    methods (Access = protected)

        function dvf = calculateScenarioDvf(this, ct, scen)
            % matRad DVF conventions (cf. matRad_doseAcc / matRad_addMovement):
            %   pull: field on the reference grid, refCube(pos) corresponds
            %         to scenCube(pos - dvf(pos))
            %   push: field on the scenario grid, scenCube(pos) corresponds
            %         to refCube(pos + dvf(pos))
            % imregdemons returns a field D with fixed(pos) ~ moving(pos + D(pos)),
            % so the pull field needs a sign flip while push can be used directly.
            switch this.dvfType
                case 'pull'
                    movingCube = ct.cubeHU{scen};
                    fixedCube = ct.cubeHU{this.refScen};
                    dvfSign = -1;
                case 'push'
                    movingCube = ct.cubeHU{this.refScen};
                    fixedCube = ct.cubeHU{scen};
                    dvfSign = 1;
            end

            registrationDvf = imregdemons(movingCube, fixedCube, this.numIterations, ...
                                          'PyramidLevels', this.pyramidLevels, ...
                                          'AccumulatedFieldSmoothing', this.accumulatedFieldSmoothing, ...
                                          'DisplayWaitbar', false);

            dvf = dvfSign * permute(registrationDvf, [4 1 2 3]);
            dvf = this.convertDvfFromVoxelUnits(dvf, ct, this.dvfUnits);
        end

        function checkDependencies(this, requireWarp)
            matRad_cfg = MatRad_Config.instance();

            if matRad_cfg.isOctave
                matRad_cfg.dispError('%s requires MATLAB with the Image Processing Toolbox (imregdemons is not available in Octave).', class(this));
            end

            if ~matRad_checkEnvImageProcessingRequirements()
                matRad_cfg.dispError('%s requires an available Image Processing Toolbox license.', class(this));
            end

            if exist('imregdemons', 'file') ~= 2
                matRad_cfg.dispError('%s requires imregdemons from the Image Processing Toolbox.', class(this));
            end

            if requireWarp && exist('imwarp', 'file') ~= 2
                matRad_cfg.dispError('Contour propagation requires imwarp from the Image Processing Toolbox.');
            end
        end

    end
end
