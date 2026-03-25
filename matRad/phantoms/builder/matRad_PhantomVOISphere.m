classdef  matRad_PhantomVOISphere < matRad_PhantomVOIVolume
    % matRad_PhantomVOISphere implements a class that helps to create spheric VOIs
    %
    % References
    %     -
    %
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
    end

    methods (Access = public)

        function obj = matRad_PhantomVOISphere(name, type, radius, varargin)
            p = inputParser;
            addParameter(p, 'objectives', {});
            addParameter(p, 'offset', [0, 0, 0]);
            addParameter(p, 'HU', 0);
            addParameter(p, 'coordType', 'voxel', @(x) numel(validatestring(x, {'voxel', 'mm'}))); % numel trick to guarantee logical cast
            parse(p, varargin{:});

            obj@matRad_PhantomVOIVolume(name, type, p); % call superclass constructor
            obj.radius = radius;
        end

        function [cst] = initializeParameters(obj, ct, cst)
            % add this VOI to the phantomBuilders cst
            ct = matRad_getWorldAxes(ct);
            cst = initializeParameters@matRad_PhantomVOIVolume(obj, cst);

            % Swaps [i j k] (x-first) <-> [j i k] (y-first / MATLAB array order)
            dimPerm = [0 1 0; 1 0 0; 0 0 1];

            % center as continuuos [j i k]
            centerPoint = (ct.cubeDim + 1) / 2;

            switch obj.coordType
                case 'voxel'
                    % Grid in [j i k]: y (rows) along dim1, x (cols) along dim2
                    [y, x, z] = ndgrid(1:ct.cubeDim(1), 1:ct.cubeDim(2), 1:ct.cubeDim(3));

                case 'mm'
                    % cubeIndex2worldCoords expects [j i k], outputs [x y z];
                    % apply dimPerm to arrive at [y x z] = [j i k] in world mm
                    centerPoint = matRad_cubeIndex2worldCoords(centerPoint, ct) * dimPerm;
                    % ct.y has nRows elements (dim1), ct.x has nCols elements (dim2)
                    [y, x, z] = ndgrid(ct.y, ct.x, ct.z);
            end

            % offset is always in [i j k]; convert to [j i k] before adding
            centerPoint = centerPoint + obj.offset * dimPerm;

            % Both modes: grid and center are in [j i k] - no extra permutation needed
            voiHelper = vecnorm([y(:) x(:) z(:)] - centerPoint, 2, 2) < obj.radius;
            voiHelper = reshape(voiHelper, ct.cubeDim);

            cst{end, 4}{1} = find(voiHelper);

        end

    end

    % Set Methods
    methods

        function set.radius(obj, value)
            validateattributes(value, {'numeric'}, {'scalar', 'positive'});
            obj.radius = value;
        end

    end
end
