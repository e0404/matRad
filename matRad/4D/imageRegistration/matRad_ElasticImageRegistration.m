classdef matRad_ElasticImageRegistration < matRad_ImageRegistration
    % matRad_ElasticImageRegistration Demons-based elastic image registration
    %
    % call
    %   obj = matRad_ElasticImageRegistration(ct, cst)
    %   obj = matRad_ElasticImageRegistration(ct, cst, refScen)
    %   obj = matRad_ElasticImageRegistration(ct, cst, refScen, metadata)
    %   obj = matRad_ElasticImageRegistration(dataStruct)
    %
    % input
    %   ct:             matRad ct struct with 4D CT scenarios
    %   cst:            matRad cst struct
    %   refScen:        reference CT scenario index (default: 1)
    %   metadata:       struct with optional registration settings
    %                   metadata.dvfType:                   'pull' or 'push' (default: 'pull')
    %                   metadata.dvfUnits:                  'voxel' or 'mm' (default: 'voxel')
    %                   metadata.numIterations:             number of registration iterations (default: 100)
    %                   metadata.pyramidLevels:             number of pyramid levels (default: 1)
    %                   metadata.accumulatedFieldSmoothing: accumulated field smoothing (default: 1)
    %
    % output
    %   obj:            elastic image registration object
    %
    % References
    %   -
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
        name = 'Elastic Registration'
    end

    properties
        ct
        cst
        refScen
        metadata
    end

    methods

        function obj = matRad_ElasticImageRegistration(ct, cst, refScen, metadata)
            obj@matRad_ImageRegistration();

            if nargin == 1 && isstruct(ct)
                obj = obj.assignPropertiesFromStruct(ct);
                return
            end

            if nargin < 2
                matRadCfg = MatRad_Config.instance();
                matRadCfg.dispError('ct and cst must be provided for elastic image registration.');
            end

            if nargin < 3 || isempty(refScen)
                refScen = 1;
            end

            if nargin < 4 || isempty(metadata)
                metadata = struct();
            end

            matRad_ElasticImageRegistration.validateCt(ct);
            refScen = matRad_ElasticImageRegistration.validateReferenceScenario(ct, refScen);
            metadata = matRad_ElasticImageRegistration.normalizeMetadata(metadata);

            obj.ct = ct;
            obj.cst = cst;
            obj.refScen = refScen;
            obj.metadata = metadata;
        end

        function ct = calcDVF(obj)
            matRadCfg = MatRad_Config.instance();
            matRad_ElasticImageRegistration.checkImageProcessingDependency(false);
            matRad_ElasticImageRegistration.validateCt(obj.ct);

            numOfCtScen = obj.ct.numOfCtScen;
            obj.ct.dvf = cell(1, numOfCtScen);
            obj.ct.dvfType = obj.metadata.dvfType;
            obj.ct.dvfUnits = obj.metadata.dvfUnits;
            obj.ct.dvfDim = obj.ct.cubeDim;
            obj.ct.refScen = obj.refScen;
            obj.ct.dvfMetadata = obj.metadata;
            obj.ct.dvfMetadata.refScen = obj.refScen;
            obj.ct.dvfMetadata.referenceCtScen = obj.refScen;

            for scen = 1:numOfCtScen
                matRadCfg.dispInfo('Registering CT scenario %d.\n', scen);
                if scen == obj.refScen
                    dvf = zeros([3 obj.ct.cubeDim]);
                else
                    dvf = obj.calculateScenarioDvf(scen);
                end
                obj.ct.dvf{scen} = dvf;
            end

            ct = obj.ct;
        end

        function [ct, cst] = propContours(obj)
            matRadCfg = MatRad_Config.instance();
            matRad_ElasticImageRegistration.checkImageProcessingDependency(true);

            if ~strcmp(obj.metadata.dvfType, 'push')
                matRadCfg.dispError('Contour propagation requires push DVFs.');
            end

            if ~isfield(obj.ct, 'dvf') || isempty(obj.ct.dvf)
                obj.ct = obj.calcDVF();
            end

            obj.ct = matRad_ElasticImageRegistration.validateDvfContainer(obj.ct, matRadCfg);
            dvfUnits = matRad_ElasticImageRegistration.getDvfUnits(obj.ct, obj.metadata);

            for structure = 1:size(obj.cst, 1)
                matRad_ElasticImageRegistration.validateReferenceContour(obj.cst, structure, obj.refScen, matRadCfg);
                refContour = obj.cst{structure, 4}{1, obj.refScen};

                if isempty(refContour)
                    matRadCfg.dispWarning('Structure %d has no reference contour. Leaving propagated contours empty.\n', structure);
                    obj.cst = matRad_ElasticImageRegistration.clearNonReferenceContours(obj.cst, structure, obj.refScen, obj.ct.numOfCtScen);
                    continue
                end

                referenceMask = false(obj.ct.cubeDim);
                referenceMask(refContour) = true;
                matRadCfg.dispInfo('Propagating contours of structure %d.\n', structure);

                for scen = 1:obj.ct.numOfCtScen
                    if scen == obj.refScen
                        continue
                    end

                    scenarioDvf = matRad_ElasticImageRegistration.convertDvfToVoxelUnits(obj.ct.dvf{scen}, obj.ct, dvfUnits, matRadCfg);
                    warpedMask = imwarp(single(referenceMask), permute(scenarioDvf, [2 3 4 1]));
                    obj.cst{structure, 4}{1, scen} = find(warpedMask > 0.5);
                end
            end

            ct = obj.ct;
            cst = obj.cst;
        end

    end

    methods (Access = private)

        function obj = assignPropertiesFromStruct(obj, inputStruct)
            names = fieldnames(inputStruct);
            for i = 1:numel(names)
                try
                    obj.(names{i}) = inputStruct.(names{i});
                catch
                end
            end
        end

        function dvf = calculateScenarioDvf(obj, scen)
            switch obj.metadata.dvfType
                case 'push'
                    movingCube = obj.ct.cubeHU{obj.refScen};
                    fixedCube = obj.ct.cubeHU{scen};
                case 'pull'
                    movingCube = obj.ct.cubeHU{scen};
                    fixedCube = obj.ct.cubeHU{obj.refScen};
            end

            registrationDvf = imregdemons(movingCube, fixedCube, obj.metadata.numIterations, ...
                                          'PyramidLevels', obj.metadata.pyramidLevels, ...
                                          'AccumulatedFieldSmoothing', obj.metadata.accumulatedFieldSmoothing);
            dvf = permute(registrationDvf, [4 1 2 3]);
            dvf = matRad_ElasticImageRegistration.convertDvfFromVoxelUnits(dvf, obj.ct, obj.metadata.dvfUnits);
        end

    end

    methods (Static, Access = private)

        function metadata = normalizeMetadata(metadata)
            matRadCfg = MatRad_Config.instance();

            if ~isstruct(metadata)
                matRadCfg.dispError('metadata must be a struct.');
            end

            metadata = matRad_ElasticImageRegistration.setDefaultMetadataField(metadata, 'dvfType', 'pull');
            metadata = matRad_ElasticImageRegistration.setDefaultMetadataField(metadata, 'dvfUnits', 'voxel');
            metadata = matRad_ElasticImageRegistration.setDefaultMetadataField(metadata, 'numIterations', 100);
            metadata = matRad_ElasticImageRegistration.setDefaultMetadataField(metadata, 'pyramidLevels', 1);
            metadata = matRad_ElasticImageRegistration.setDefaultMetadataField(metadata, 'accumulatedFieldSmoothing', 1);

            metadata.dvfType = char(metadata.dvfType);
            metadata.dvfUnits = char(metadata.dvfUnits);

            if ~any(strcmp(metadata.dvfType, {'pull', 'push'}))
                matRadCfg.dispError('metadata.dvfType must be ''pull'' or ''push''.');
            end

            if ~any(strcmp(metadata.dvfUnits, {'voxel', 'mm'}))
                matRadCfg.dispError('metadata.dvfUnits must be ''voxel'' or ''mm''.');
            end

            matRad_ElasticImageRegistration.validatePositiveInteger(metadata.numIterations, 'metadata.numIterations', matRadCfg);
            matRad_ElasticImageRegistration.validatePositiveInteger(metadata.pyramidLevels, 'metadata.pyramidLevels', matRadCfg);
            matRad_ElasticImageRegistration.validateNonnegativeScalar(metadata.accumulatedFieldSmoothing, ...
                                                                      'metadata.accumulatedFieldSmoothing', matRadCfg);
        end

        function metadata = setDefaultMetadataField(metadata, fieldName, defaultValue)
            if ~isfield(metadata, fieldName) || isempty(metadata.(fieldName))
                metadata.(fieldName) = defaultValue;
            end
        end

        function checkImageProcessingDependency(requireWarp)
            matRadCfg = MatRad_Config.instance();
            env = matRad_getEnvironment();

            if strcmp(env, 'OCTAVE')
                matRadCfg.dispError('matRad_ElasticImageRegistration requires MATLAB with the Image Processing Toolbox.');
            end

            if exist('imregdemons', 'file') ~= 2
                matRadCfg.dispError('matRad_ElasticImageRegistration requires imregdemons from the Image Processing Toolbox.');
            end

            if requireWarp && exist('imwarp', 'file') ~= 2
                matRadCfg.dispError('Contour propagation requires imwarp from the Image Processing Toolbox.');
            end

            if exist('license', 'builtin') == 5 || exist('license', 'file') == 2
                hasImageToolboxLicense = true;
                try
                    hasImageToolboxLicense = license('test', 'image_toolbox');
                catch
                end

                if ~hasImageToolboxLicense
                    matRadCfg.dispError('matRad_ElasticImageRegistration requires an available Image Processing Toolbox license.');
                end
            end
        end

        function validateCt(ct)
            matRadCfg = MatRad_Config.instance();

            if ~isstruct(ct)
                matRadCfg.dispError('ct must be a struct.');
            end

            if ~isfield(ct, 'numOfCtScen') || isempty(ct.numOfCtScen)
                matRadCfg.dispError('ct.numOfCtScen must be defined.');
            end

            matRad_ElasticImageRegistration.validatePositiveInteger(ct.numOfCtScen, 'ct.numOfCtScen', matRadCfg);

            if ~isfield(ct, 'cubeDim') || numel(ct.cubeDim) ~= 3
                matRadCfg.dispError('ct.cubeDim must be a three-element vector.');
            end

            if ~isfield(ct, 'cubeHU') || ~iscell(ct.cubeHU) || numel(ct.cubeHU) < ct.numOfCtScen
                matRadCfg.dispError('ct.cubeHU must contain all CT scenarios.');
            end

            for scen = 1:ct.numOfCtScen
                if ~isequal(size(ct.cubeHU{scen}), ct.cubeDim)
                    matRadCfg.dispError('ct.cubeHU{%d} size must match ct.cubeDim.', scen);
                end
            end
        end

        function refScen = validateReferenceScenario(ct, refScen)
            matRadCfg = MatRad_Config.instance();
            matRad_ElasticImageRegistration.validatePositiveInteger(refScen, 'refScen', matRadCfg);
            refScen = double(refScen);

            if refScen > ct.numOfCtScen
                matRadCfg.dispError('refScen (%d) exceeds ct.numOfCtScen (%d).', refScen, ct.numOfCtScen);
            end
        end

        function validatePositiveInteger(value, name, matRadCfg)
            if ~(isnumeric(value) && isscalar(value) && isfinite(value) && round(value) == value && value >= 1)
                matRadCfg.dispError('%s must be a positive integer scalar.', name);
            end
        end

        function validateNonnegativeScalar(value, name, matRadCfg)
            if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0)
                matRadCfg.dispError('%s must be a nonnegative scalar.', name);
            end
        end

        function ct = validateDvfContainer(ct, matRadCfg)
            if ~isfield(ct, 'dvf') || ~iscell(ct.dvf) || numel(ct.dvf) < ct.numOfCtScen
                matRadCfg.dispError('ct.dvf must contain all CT scenario deformation vector fields.');
            end

            for scen = 1:ct.numOfCtScen
                matRad_ElasticImageRegistration.validateDvf(ct.dvf{scen}, ct, scen, matRadCfg);
            end
        end

        function validateDvf(dvf, ct, scen, matRadCfg)
            if ndims(dvf) ~= 4 || size(dvf, 1) ~= 3 || ~isequal(size(dvf), [3 ct.cubeDim])
                matRadCfg.dispError('ct.dvf{%d} must have size [3 ct.cubeDim].', scen);
            end
        end

        function validateReferenceContour(cst, structure, refScen, matRadCfg)
            if size(cst, 2) < 4 || isempty(cst{structure, 4}) || numel(cst{structure, 4}) < refScen
                matRadCfg.dispError('Structure %d does not contain contours for reference CT scenario %d.', structure, refScen);
            end
        end

        function cst = clearNonReferenceContours(cst, structure, refScen, numOfCtScen)
            for scen = 1:numOfCtScen
                if scen ~= refScen
                    cst{structure, 4}{1, scen} = [];
                end
            end
        end

        function dvfUnits = getDvfUnits(ct, metadata)
            dvfUnits = metadata.dvfUnits;
            if isfield(ct, 'dvfUnits') && ~isempty(ct.dvfUnits)
                dvfUnits = char(ct.dvfUnits);
            elseif isfield(ct, 'dvfMetadata') && isfield(ct.dvfMetadata, 'dvfUnits') && ~isempty(ct.dvfMetadata.dvfUnits)
                dvfUnits = char(ct.dvfMetadata.dvfUnits);
            end
        end

        function dvf = convertDvfFromVoxelUnits(dvf, ct, dvfUnits)
            if strcmp(dvfUnits, 'mm')
                matRad_ElasticImageRegistration.validateResolution(ct, MatRad_Config.instance());
                dvf(1, :, :, :) = dvf(1, :, :, :) .* ct.resolution.x;
                dvf(2, :, :, :) = dvf(2, :, :, :) .* ct.resolution.y;
                dvf(3, :, :, :) = dvf(3, :, :, :) .* ct.resolution.z;
            end
        end

        function dvf = convertDvfToVoxelUnits(dvf, ct, dvfUnits, matRadCfg)
            switch dvfUnits
                case 'voxel'
                    return
                case 'mm'
                    matRad_ElasticImageRegistration.validateResolution(ct, matRadCfg);
                    dvf(1, :, :, :) = dvf(1, :, :, :) ./ ct.resolution.x;
                    dvf(2, :, :, :) = dvf(2, :, :, :) ./ ct.resolution.y;
                    dvf(3, :, :, :) = dvf(3, :, :, :) ./ ct.resolution.z;
                otherwise
                    matRadCfg.dispError('Unsupported DVF units: %s.', dvfUnits);
            end
        end

        function validateResolution(ct, matRadCfg)
            if ~isfield(ct, 'resolution') || ~isfield(ct.resolution, 'x') || ~isfield(ct.resolution, 'y') || ~isfield(ct.resolution, 'z')
                matRadCfg.dispError('ct.resolution.x, ct.resolution.y and ct.resolution.z are required for DVFs in mm.');
            end
        end

    end
end
