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
        boxDimensionsType  % either voxel or mm
    end

    methods (Access = public)

        function obj = matRad_PhantomVOIBox(name, type, boxDimensions, varargin)
            p = inputParser;
            addParameter(p, 'objectives', {});
            addParameter(p, 'offset', [0, 0, 0]);
            addParameter(p, 'HU', 0);
            addParameter(p, 'boxDimensionsType', 'voxel', @(x) numel(validatestring(x, {'voxel', 'mm'})));
            parse(p, varargin{:});

            obj@matRad_PhantomVOIVolume(name, type, p); % call superclass constructor
            obj.boxDimensions = boxDimensions;
            obj.boxDimensionsType = p.Results.boxDimensionsType;
        end

        function [cst] = initializeParameters(obj, ct, cst)
            % add this objective to the phantomBuilders cst
            ct = matRad_getWorldAxes(ct);
            cst = initializeParameters@matRad_PhantomVOIVolume(obj, cst);

            dimPerm = [0 1 0; 1 0 0; 0 0 1];

            % center = round(ct.cubeDim/2);
            % from here we work with continuous voxel indices
            center = ct.cubeDim / 2;
            offsets = obj.offset;
            dims = obj.boxDimensions;
            centerPoint = center * dimPerm + offsets;

            switch obj.boxDimensionsType
                case 'mm' % obtain points and limits in world coordinates
                    centerPoint = matRad_cubeIndex2worldCoords(centerPoint, ct);
                    halfRes =  [ct.resolution.x ct.resolution.y ct.resolution.z] / 2;
                    ctMin = [min(ct.x) min(ct.y) min(ct.z)] - halfRes;
                    ctMax = [max(ct.x) max(ct.y) max(ct.z)] + halfRes;

                    [x, y, z] = ndgrid(ct.x, ct.y, ct.z);
                case 'voxel' % obtain limits in voxel indices
                    ctMin = [1 1 1];
                    ctMax = ct.cubeDim * dimPerm;
                    [y, x, z] = ndgrid(1:ct.cubeDim(1), 1:ct.cubeDim(2), 1:ct.cubeDim(3));
            end
            coords = [x(:) y(:) z(:)];

            maxPoints = min(centerPoint + dims / 2, ctMax);
            minPoints = max(centerPoint - dims / 2, ctMin);

            voiHelper = all(coords >= minPoints & coords <= maxPoints, 2);
            voiHelper = reshape(voiHelper, ct.cubeDim);

            cst{end, 4}{1} = find(voiHelper);
        end

    end

    % Set Methods
    methods

        function set.boxDimensionsType(obj, dimType)
            obj.boxDimensionsType = validatestring(dimType, {'voxel', 'mm'});
        end

        function set.boxDimensions(obj, dims)
            validateattributes(dims, {'numeric'}, {'vector', 'numel', 3, 'positive'});
            obj.boxDimensions = dims;
        end

    end
end
