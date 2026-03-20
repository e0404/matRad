classdef matRad_PhantomVOIBox < matRad_PhantomVOIVolume
    % matRad_PhantomVOIBox implements a class that helps to create box VOIs
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
    properties % additional property of cubic objects
        boxDimensions
    end

    methods (Access = public)

        function obj = matRad_PhantomVOIBox(name, type, boxDimensions, varargin)
            p = inputParser;
            addParameter(p, 'objectives', {});
            addParameter(p, 'offset', [0, 0, 0]);
            addParameter(p, 'HU', 0);
            addParameter(p, 'coordType', 'voxel', @(x) numel(validatestring(x, {'voxel', 'mm'})));
            parse(p, varargin{:});

            obj@matRad_PhantomVOIVolume(name, type, p); % call superclass constructor
            obj.boxDimensions = boxDimensions;
        end

        function [cst] = initializeParameters(obj, ct, cst)
            % add this objective to the phantomBuilders cst
            ct = matRad_getWorldAxes(ct);
            cst = initializeParameters@matRad_PhantomVOIVolume(obj, cst);

            % Swaps [i j k] (x-first) <-> [j i k] (y-first / MATLAB array order)
            dimPerm = [0 1 0; 1 0 0; 0 0 1];

            % Center in [j i k] (cubeDim is already in MATLAB array order)
            centerPoint = (ct.cubeDim + 1) / 2;

            switch obj.coordType
                case 'voxel'
                    ctMin = [1 1 1];
                    ctMax = ct.cubeDim;  % [j i k]
                    [y, x, z] = ndgrid(1:ct.cubeDim(1), 1:ct.cubeDim(2), 1:ct.cubeDim(3));

                case 'mm'
                    % cubeIndex2worldCoords expects [i j k], outputs [x y z];
                    % * dimPerm converts to [y x z] = [j i k] in world mm
                    centerPoint = matRad_cubeIndex2worldCoords(centerPoint, ct) * dimPerm;
                    halfRes = [ct.resolution.y ct.resolution.x ct.resolution.z] / 2;
                    ctMin = [min(ct.y) min(ct.x) min(ct.z)] - halfRes;
                    ctMax = [max(ct.y) max(ct.x) max(ct.z)] + halfRes;
                    % ct.y has nRows elements (dim1), ct.x has nCols elements (dim2)
                    [y, x, z] = ndgrid(ct.y, ct.x, ct.z);
            end

            % offset and boxDimensions are in [i j k]; convert to [j i k]
            centerPoint = centerPoint + obj.offset * dimPerm;
            dims = obj.boxDimensions * dimPerm;

            coords = [y(:) x(:) z(:)];  % [j i k]

            maxPoints = min(centerPoint + dims / 2, ctMax);
            minPoints = max(centerPoint - dims / 2, ctMin);

            voiHelper = all(coords >= minPoints & coords <= maxPoints, 2);
            voiHelper = reshape(voiHelper, ct.cubeDim);

            cst{end, 4}{1} = find(voiHelper);
        end

    end

    % Set Methods
    methods

        function set.boxDimensions(obj, dims)
            validateattributes(dims, {'numeric'}, {'vector', 'numel', 3, 'positive'});
            obj.boxDimensions = dims;
        end

    end
end
