classdef (Abstract) matRad_ImageRegistrationBase < handle
    % matRad_ImageRegistrationBase Base class for 4D image registration
    %   Defines the common interface for image registration classes that
    %   compute deformation vector fields (DVFs) between CT scenarios and
    %   propagate contours between them.
    %
    %   DVF conventions (compatible with matRad_doseAcc / matRad_addMovement):
    %
    %   - ct.dvf{scen} has size [3 ct.cubeDim] with component order x,y,z,
    %     where x refers to the second cube dimension (columns), y to the
    %     first (rows) and z to the third (slices).
    %   - 'pull' DVFs are defined on the reference scenario grid such that
    %     refCube(pos) corresponds to scenCube(pos - dvf(pos)).
    %   - 'push' DVFs are defined on the scenario grid such that
    %     scenCube(pos) corresponds to refCube(pos + dvf(pos)).
    %   - dvfUnits is 'mm' by default, as expected by matRad_doseAcc.
    %   - Registration settings used to create the DVFs are stored in
    %     ct.dvfMetadata (fields dvfType, dvfUnits, refScen).
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

    properties (Abstract, Constant)
        name % user readable name of the registration algorithm
    end

    properties
        refScen  = 1       % reference CT scenario index
        dvfType  = 'pull'  % 'pull' or 'push' deformation vector fields
        dvfUnits = 'mm'    % 'mm' or 'voxel' - matRad's 4D workflows (matRad_doseAcc) expect [mm]
    end

    methods

        function this = matRad_ImageRegistrationBase(propStruct)
            if nargin > 0 && ~isempty(propStruct)
                this.assignProperties(propStruct);
            end
        end

        function set.refScen(this, value)
            matRad_cfg = MatRad_Config.instance();
            if ~(isnumeric(value) && isscalar(value) && isfinite(value) && round(value) == value && value >= 1)
                matRad_cfg.dispError('refScen must be a positive integer scalar.');
            end
            this.refScen = double(value);
        end

        function set.dvfType(this, value)
            matRad_cfg = MatRad_Config.instance();
            if ~((ischar(value) || isstring(value)) && any(strcmp(char(value), {'pull', 'push'})))
                matRad_cfg.dispError('dvfType must be ''pull'' or ''push''.');
            end
            this.dvfType = char(value);
        end

        function set.dvfUnits(this, value)
            matRad_cfg = MatRad_Config.instance();
            if ~((ischar(value) || isstring(value)) && any(strcmp(char(value), {'mm', 'voxel'})))
                matRad_cfg.dispError('dvfUnits must be ''mm'' or ''voxel''.');
            end
            this.dvfUnits = char(value);
        end

        function s = struct(this)
            % Serializes the registration settings (not any patient data)
            s = struct('className', class(this));
            props = properties(this);
            for i = 1:numel(props)
                if ~strcmp(props{i}, 'name')
                    s.(props{i}) = this.(props{i});
                end
            end
        end

    end

    % Would be abstract methods, but Octave does not support abstract method
    % declarations, so we define them with an error (cf. matRad_DoseEngineBase)
    methods

        function ct = calcDVF(this, ct)
            % Computes deformation vector fields between the reference CT
            % scenario and all other scenarios and stores them in
            % ct.dvf / ct.dvfMetadata
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('calcDVF needs to be implemented by the used image registration subclass!');
        end

        function [ct, cst] = propContours(this, ct, cst)
            % Propagates the contours of the reference CT scenario to all
            % other scenarios using push deformation vector fields
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('propContours needs to be implemented by the used image registration subclass!');
        end

    end

    methods (Access = protected)

        function assignProperties(this, propStruct)
            % Assigns properties from a (serialized) struct. Unknown fields
            % are reported as warnings, invalid values throw
            matRad_cfg = MatRad_Config.instance();

            if ~isstruct(propStruct)
                matRad_cfg.dispError('Registration properties must be provided as struct.');
            end

            fields = fieldnames(propStruct);
            settableProps = properties(this);
            for i = 1:numel(fields)
                if any(strcmp(fields{i}, {'className', 'name'}))
                    continue
                end
                if ~any(strcmp(fields{i}, settableProps))
                    matRad_cfg.dispWarning('Ignoring unknown registration property ''%s'' for %s.', fields{i}, class(this));
                    continue
                end
                this.(fields{i}) = propStruct.(fields{i});
            end
        end

        function validateCt(this, ct)
            matRad_cfg = MatRad_Config.instance();

            if ~isstruct(ct)
                matRad_cfg.dispError('ct must be a struct.');
            end

            if ~isfield(ct, 'numOfCtScen') || ...
                ~(isnumeric(ct.numOfCtScen) && isscalar(ct.numOfCtScen) && ct.numOfCtScen >= 1 && round(ct.numOfCtScen) == ct.numOfCtScen)
                matRad_cfg.dispError('ct.numOfCtScen must be a positive integer scalar.');
            end

            if ~isfield(ct, 'cubeDim') || numel(ct.cubeDim) ~= 3
                matRad_cfg.dispError('ct.cubeDim must be a three-element vector.');
            end

            if ~isfield(ct, 'cubeHU') || ~iscell(ct.cubeHU) || numel(ct.cubeHU) < ct.numOfCtScen
                matRad_cfg.dispError('ct.cubeHU must contain all CT scenarios.');
            end

            for scen = 1:ct.numOfCtScen
                if ~isequal(size(ct.cubeHU{scen}), reshape(ct.cubeDim, 1, []))
                    matRad_cfg.dispError('ct.cubeHU{%d} size must match ct.cubeDim.', scen);
                end
            end

            if this.refScen > ct.numOfCtScen
                matRad_cfg.dispError('refScen (%d) exceeds ct.numOfCtScen (%d).', this.refScen, ct.numOfCtScen);
            end
        end

        function validateResolution(this, ct)
            matRad_cfg = MatRad_Config.instance();
            if ~isfield(ct, 'resolution') || ~all(isfield(ct.resolution, {'x', 'y', 'z'}))
                matRad_cfg.dispError('ct.resolution.x/y/z are required for DVFs in mm.');
            end
        end

        function validateDvfContainer(this, ct)
            matRad_cfg = MatRad_Config.instance();

            if ~isfield(ct, 'dvf') || ~iscell(ct.dvf) || numel(ct.dvf) < ct.numOfCtScen
                matRad_cfg.dispError('ct.dvf must contain all CT scenario deformation vector fields.');
            end

            for scen = 1:ct.numOfCtScen
                if ndims(ct.dvf{scen}) ~= 4 || ~isequal(size(ct.dvf{scen}), [3 reshape(ct.cubeDim, 1, [])])
                    matRad_cfg.dispError('ct.dvf{%d} must have size [3 ct.cubeDim].', scen);
                end
            end
        end

        function metadata = getDvfMetadata(this, ct)
            % Effective DVF metadata: taken from ct.dvfMetadata if the DVFs
            % travel with the ct, completed with this object's settings
            metadata = struct('dvfType', this.dvfType, ...
                              'dvfUnits', this.dvfUnits, ...
                              'refScen', this.refScen);

            if isfield(ct, 'dvfMetadata')
                fields = fieldnames(ct.dvfMetadata);
                for i = 1:numel(fields)
                    if ~isempty(ct.dvfMetadata.(fields{i}))
                        metadata.(fields{i}) = ct.dvfMetadata.(fields{i});
                    end
                end
            end
        end

        function dvf = convertDvfToVoxelUnits(this, dvf, ct, dvfUnits)
            matRad_cfg = MatRad_Config.instance();
            switch dvfUnits
                case 'voxel'
                    % nothing to do
                case 'mm'
                    this.validateResolution(ct);
                    dvf(1, :, :, :) = dvf(1, :, :, :) ./ ct.resolution.x;
                    dvf(2, :, :, :) = dvf(2, :, :, :) ./ ct.resolution.y;
                    dvf(3, :, :, :) = dvf(3, :, :, :) ./ ct.resolution.z;
                otherwise
                    matRad_cfg.dispError('Unsupported DVF units: %s.', dvfUnits);
            end
        end

        function dvf = convertDvfFromVoxelUnits(this, dvf, ct, dvfUnits)
            matRad_cfg = MatRad_Config.instance();
            switch dvfUnits
                case 'voxel'
                    % nothing to do
                case 'mm'
                    this.validateResolution(ct);
                    dvf(1, :, :, :) = dvf(1, :, :, :) .* ct.resolution.x;
                    dvf(2, :, :, :) = dvf(2, :, :, :) .* ct.resolution.y;
                    dvf(3, :, :, :) = dvf(3, :, :, :) .* ct.resolution.z;
                otherwise
                    matRad_cfg.dispError('Unsupported DVF units: %s.', dvfUnits);
            end
        end

    end
end
