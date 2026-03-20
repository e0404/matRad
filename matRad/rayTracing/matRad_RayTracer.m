classdef matRad_RayTracer < handle
    % matRad_RayTracer Base class for ray tracing through a patient CT volume
    %   This class provides the infrastructure for ray tracing in beam's-eye-view
    %   (BEV) coordinates to compute radiological depths along rays from a
    %   radiation source through a patient CT grid.
    %
    %   Subclasses must implement the abstract method:
    %       traceRay()   - trace a single ray and return path segment lengths,
    %                      density values, and intersected voxel indices
    %   Subclasses can optionally override:
    %       traceRays()  - trace multiple rays in a vectorized or otherwise accelerated manner
    %                      (default: loop over traceRay)
    %
    %   Usage example:
    %       tracer = matRad_SomeRayTracerSubclass(cubes, grid);
    %       [radDepthsV, radDepthCube] = tracer.traceCube(stfElement, V);
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

    properties
        lateralCutOff
    end

    properties (Access = protected)
        cubes
        planes
        grid
        numPlanes
    end

    properties (Hidden)
        useOldCandidateRayMatrixCalculation = false
    end

    methods

        function this = matRad_RayTracer(cubes, grid)
            % matRad_RayTracer Construct a ray tracer instance
            %   cubes - cell array of density/material cubes to trace through
            %   grid  - CT grid struct with fields x, y, z and resolution

            this.cubes = cubes;
            this.grid  = grid;

            this.initializeGeometry();

            matRad_cfg = MatRad_Config.instance();
            this.lateralCutOff = matRad_cfg.defaults.propDoseCalc.geometricLateralCutOff;
        end

        function [alphas, l, rhoTmp, d12, ix] = traceRays(this, ...
                                                          isocenter, ...
                                                          sourcePoints, ...
                                                          targetPoints)

            % Default trivial implementation based on traceRay
            nRays = size(targetPoints, 1);
            nSources = size(sourcePoints, 1);

            if nSources ~= nRays && nSources ~= 1
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Number of source points (%d) needs to be one or equal to number of target points (%d)!', nSources, nRays);
            elseif nSources == 1
                sourcePoints = repmat(sourcePoints, nRays, 1);
                nSources = nRays;
            end

            nCubes = numel(this.cubes);

            rhoTmp = cell(nRays, nCubes);
            alphas = cell(nRays, 1);
            l = cell(nRays, 1);
            ix = cell(nRays, 1);
            d12 = NaN(nRays, 1);
            for r = 1:nRays
                [alphas{r}, l{r}, rhoTmp(r, :), d12(r), ix{r}] = this.traceRay(isocenter, sourcePoints(r, :), targetPoints(r, :));
            end

            % pad with NaN values
            maxnumalphas = max(cellfun(@numel, alphas));
            maxnumix = max(cellfun(@numel, ix));

            nanpad = @(x, maxval) [x(1:end), NaN(1, maxnumalphas - length(x))];
            alphas = cellfun(nanpad, alphas, 'UniformOutput', false);
            nanpad = @(x) [x(1:end), NaN(1, maxnumix - length(x))];

            l = cellfun(nanpad, l, 'UniformOutput', false);
            ix = cellfun(nanpad, ix, 'UniformOutput', false);
            rhoTmp = cellfun(nanpad, rhoTmp, 'UniformOutput', false);

            % now make matrices
            alphas = cell2mat(alphas');
            l = cell2mat(l');
            ix = cell2mat(ix');

            rho = cell(1, nCubes);
            for c = flip(1:nCubes)
                rho{c} = cell2mat(rhoTmp(:, c));
                rhoTmp(:, c) = [];
            end

        end

        function [alphas, l, rho, d12, ix] = traceRay(this, ...
                                                      isocenter, ...
                                                      sourcePoint, ...
                                                      targetPoint)
            error('Needs to be implemented!');
        end

        function [radDepthsV, radDepthCube] = traceCube(this, stfElement, voxelIndices, rotCoordsV)
            matRad_cfg = MatRad_Config.instance();

            if ~isstruct(stfElement) || numel(stfElement) ~= 1
                matRad_cfg.dispError('The RayTracer does not accept stf struct arrays and only operates on a single field!');
            end

            cubeDim = size(this.cubes{1});
            nCubeVox = numel(this.cubes{1});

            % If no subset of voxels is specified, take all of them
            if nargin < 4
                voxelIndices = transpose(1:nCubeVox);
                if matRad_cfg.enableGPU
                    voxelIndices = gpuArray(int32(voxelIndices));
                end
            end

            % if we don't provide rotated patient coordinates we can compute
            % them here on our own
            if nargin < 5
                coordsV = single(matRad_cubeIndex2worldCoords(voxelIndices, this.grid));

                % Get Rotation Matrix
                % Do not transpose matrix since we usage of row vectors &
                % transformation of the coordinate system need double transpose

                rotMat_system_T = matRad_getRotationMatrix(stfElement.gantryAngle, stfElement.couchAngle);

                % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
                rotCoordsV         = (coordsV - stfElement.isoCenter) * rotMat_system_T;

                % translate relative to source
                rotCoordsV = rotCoordsV - stfElement.sourcePoint_bev;
            end

            % set up ray matrix direct behind last voxel
            rayMx_bev_y = max(rotCoordsV(:, 2)) + max([this.grid.resolution.x this.grid.resolution.y this.grid.resolution.z]);
            rayMx_bev_y = rayMx_bev_y + stfElement.sourcePoint_bev(2);
            rayMatrixScale = 1 + rayMx_bev_y / stfElement.SAD;

            % set up list with bev coordinates for calculation of radiological depth
            coords = zeros(nCubeVox, 3, class(rotCoordsV));
            coords(voxelIndices, :) = rotCoordsV;

            referencePositionsBEV = rayMatrixScale * vertcat(stfElement.ray.rayPos_bev);

            % calculate spacing of rays on ray matrix
            rayMxSpacing = cast(1 / sqrt(2) * ...
                                min([this.grid.resolution.x this.grid.resolution.y this.grid.resolution.z]), matRad_underlyingTypeCompat(coords));
            spacingRange = rayMxSpacing * (floor(-500 / rayMxSpacing):ceil(500 / rayMxSpacing));

            t_candidate = tic;
            candidateRayMx = this.getCandidateRayMatrix(spacingRange, referencePositionsBEV);
            t_candidate = toc(t_candidate);

            [rayIdxZ, rayIdxX] = find(candidateRayMx);
            rayMx_bev = [spacingRange(rayIdxX)' rayMx_bev_y * ones(sum(candidateRayMx(:)), 1) spacingRange(rayIdxZ)'];

            % Rotation matrix. Transposed because of row vectors
            rotMat_vectors_T = transpose(matRad_getRotationMatrix(stfElement.gantryAngle, stfElement.couchAngle));

            % rotate ray matrix from bev to world coordinates
            rayMx_world = rayMx_bev * rotMat_vectors_T;

            % criterium for ray selection
            raySelection = rayMxSpacing / 2;

            % Trace all selected rays
            [~, l, rho, ~, ixHitVoxel] = this.traceRays(stfElement.isoCenter, stfElement.sourcePoint, rayMx_world);

            % find voxels for which we should remember this tracing because this is
            % the closest ray by projecting the voxel coordinates to the
            % intersection points with the ray matrix and checking if the distance
            % in x and z direction is smaller than the resolution of the ray matrix
            scale_factor = NaN(size(ixHitVoxel));
            valid_ix = ~isnan(ixHitVoxel);
            scale_factor(valid_ix) = (rayMx_bev_y - stfElement.sourcePoint_bev(2)) ./ ...
                   coords(ixHitVoxel(valid_ix), 2);

            x_dist = NaN(size(ixHitVoxel));
            z_dist = NaN(size(ixHitVoxel));

            x_dist(valid_ix) = coords(ixHitVoxel(valid_ix), 1) .* scale_factor(valid_ix);
            x_dist = x_dist - rayMx_bev(:, 1);

            z_dist(valid_ix) = coords(ixHitVoxel(valid_ix), 3) .* scale_factor(valid_ix);
            z_dist = z_dist - rayMx_bev(:, 3);

            % Find indices
            ixRememberFromCurrTracing = x_dist > -raySelection & x_dist <= raySelection & ...
                z_dist > -raySelection & z_dist <= raySelection;

            % set up rad depth cube for results
            radDepthCube = repmat({NaN(cubeDim, class(l))}, numel(this.cubes));
            radDepthsV = cell(size(radDepthCube));

            for j = 1:numel(this.cubes)
                rayDistances = l .* rho{j};
                rayWepl = cumsum(rayDistances, 2) - rayDistances / 2;

                radDepthCube{j}(ixHitVoxel(ixRememberFromCurrTracing)) = rayWepl(ixRememberFromCurrTracing);

                radDepthsV{j} = radDepthCube{j}(voxelIndices);
            end
        end

    end

    methods (Access = protected)

        function candidateRayMx = getCandidateRayMatrix(this, spacingRange, refPosBEV)
            % define candidate ray matrix covering 1000x1000mm^2
            numOfCandidateRays = numel(spacingRange);
            candidateRayMx     = zeros(numOfCandidateRays, 'logical');

            rSq = this.lateralCutOff^2;

            if this.useOldCandidateRayMatrixCalculation
                % old implementation (based on meshgrid), which is slightly
                % faster for GPUs
                [candidateRayCoordsX, candidateRayCoordsZ] = meshgrid(spacingRange);

                % check which rays should be used
                for i = 1:size(refPosBEV, 1)
                    ixCandidates = (candidateRayCoordsX(:) - refPosBEV(i, 1)).^2 + ...
                        (candidateRayCoordsZ(:) - refPosBEV(i, 3)).^2 <= rSq;

                    candidateRayMx(ixCandidates) = 1;
                end
            end

            % check which rays should be used
            for i = 1:size(refPosBEV, 1)
                % Old implementation (based on meshgrid)
                % xCandidateCoords = (candidateRayCoords_X(:) - refPosBEV(i,1)).^2;
                % zCandidateCoords = (candidateRayCoords_Z(:) - refPosBEV(i,3)).^2;

                % ixCandidates = (xCandidateCoords + zCandidateCoords) <= rSq;
                % candidateRayMx(ixCandidates) = true;

                % simple boolean update
                xCandidateCoords = (spacingRange - refPosBEV(i, 1)).^2;
                zCandidateCoords = (spacingRange - refPosBEV(i, 3)).^2;

                zOk = zCandidateCoords <= rSq;
                if ~any(zOk)
                    continue
                end

                z0 = find(zOk, 1, 'first');
                z1 = find(zOk, 1, 'last');
                candidateRayMx(z0:z1, :) = candidateRayMx(z0:z1, :) | zCandidateCoords(z0:z1)' + xCandidateCoords <= rSq;
            end
        end

        function initializeGeometry(this)
            this.grid = matRad_getWorldAxes(this.grid);

            this.planes.x = [this.grid.x - this.grid.resolution.x / 2, this.grid.x(end) + this.grid.resolution.x / 2];
            this.planes.y = [this.grid.y - this.grid.resolution.y / 2, this.grid.y(end) + this.grid.resolution.y / 2];
            this.planes.z = [this.grid.z - this.grid.resolution.z / 2, this.grid.z(end) + this.grid.resolution.z / 2];

            this.numPlanes = arrayfun(@(x) numel(this.planes.(x)), 'xyz');
        end

    end
end
