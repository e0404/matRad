classdef  matRad_PhantomVOISphere < matRad_PhantomVOIVolume
    % matRad_PhantomVOISphere implements a class that helps to create spheric VOIs
    %
    % References
    %     -
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2022 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        radius
        radiusType  % either voxel or mm
    end

    methods (Access = public)

        function obj = matRad_PhantomVOISphere(name, type, radius, varargin)
            p = inputParser;
            addParameter(p, 'objectives', {});
            addParameter(p, 'offset', [0, 0, 0]);
            addParameter(p, 'HU', 0);
            addParameter(p, 'radiusType', 'voxel', @(x) numel(validatestring(x, {'voxel', 'mm'}))); % numel trick to guarantee logical cast
            parse(p, varargin{:});

            obj@matRad_PhantomVOIVolume(name, type, p); % call superclass constructor
            obj.radius = radius;
            obj.radiusType = p.Results.radiusType;
        end

        function [cst] = initializeParameters(obj, ct, cst)
            % add this VOI to the phantomBuilders cst
            ct = matRad_getWorldAxes(ct);
            cst = initializeParameters@matRad_PhantomVOIVolume(obj, cst);
            center = ct.cubeDim / 2;
            offset = obj.offset;
            centerPoint = center - offset;
            switch obj.radiusType
                case 'voxel'
                    acceptFunc = @(currVoxel) vecnorm(currVoxel - centerPoint, 2, 2) < obj.radius;
                case 'mm'
                    centerCoord = matRad_cubeIndex2worldCoords(center, ct);
                    offsetCoord = offset .* [ct.resolution.x ct.resolution.y ct.resolution.z];
                    centerPointCoord = centerCoord - offsetCoord;
                    acceptFunc = @(currVoxel) vecnorm(matRad_cubeIndex2worldCoords(currVoxel, ct) - centerPointCoord, 2, 2) < obj.radius;
            end

            [y, x, z] = ndgrid(1:ct.cubeDim(1), 1:ct.cubeDim(2), 1:ct.cubeDim(3));
            voiHelper = acceptFunc([y(:) x(:) z(:)]);
            voiHelper = reshape(voiHelper, ct.cubeDim);

            cst{end, 4}{1} = find(voiHelper);

        end

    end
    % Set Methods
    methods

        function set.radiusType(obj, rType)
            obj.radiusType = validatestring(rType, {'voxel', 'mm'});
        end

        function set.radius(obj, value)
            validateattributes(value, {'numeric'}, {'scalar', 'positive'});
            obj.radius = value;
        end

    end
end
