classdef matRad_RayTracerSiddon < matRad_RayTracer
    % matRad_RayTracerSiddon Ray tracer using Siddon's algorithm
    %   Implements ray tracing through a CT volume using Siddon's algorithm,
    %   which analytically computes the exact intersection lengths of a ray
    %   with each voxel by tracking parametric alpha values at which the ray
    %   crosses each set of parallel CT planes (x, y, z).
    %
    %   This class overrides traceRays() from matRad_RayTracer with a
    %   vectorized implementation that handles multiple rays simultaneously,
    %   significantly improving performance over the default loop.
    %
    %   References:
    %       Siddon RL (1985), "Fast calculation of the exact radiological
    %       path for a three-dimensional CT array", Med. Phys. 12(2):252-5.
    %
    %   Usage example:
    %       tracer = matRad_RayTracerSiddon(cubes, grid);
    %       [alphas, l, rho, d12, ix] = tracer.traceRays(isocenter, sourcePoints, targetPoints);
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

    methods

        function this = matRad_RayTracerSiddon(cubes, grid)
            this = this@matRad_RayTracer(cubes, grid);
        end

        function [alphas, l, rho, d12, ix] = traceRay(this, ...
                                                      isocenter, ...
                                                      sourcePoint, ...
                                                      targetPoint)

            nRays = size(targetPoint, 1);
            nSources = size(sourcePoint, 1);

            if nRays ~= 1 || nSources ~= 1
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError(['Number of target Points and source points needs to be equal to one! ' ...
                                      'If you want to trace multiple rays at once, use traceRays instead!']);
            end

            [alphas, l, rho, d12, ix] = this.traceRays(isocenter, sourcePoint, targetPoint);
        end

        function [alphas, l, rho, d12, ix] = traceRays(this, ...
                                                       isocenter, ...
                                                       sourcePoint, ...
                                                       targetPoint)

            nRays = size(targetPoint, 1);
            nSources = size(sourcePoint, 1);

            if nSources ~= nRays && nSources ~= 1
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Number of source points (%d) needs to be one or equal to number of target points (%d)!', nSources, nRays);
            elseif nSources == 1
                sourcePoint = repmat(sourcePoint, nRays, 1);
                % nSources = nRays;
            end

            rayVec = targetPoint - sourcePoint;
            sourcePoint = sourcePoint + isocenter;

            % eq 7 & 8
            % Calculate relative distances (alphas) at which intersections
            % occur
            alphas = this.computeAllAlphas(sourcePoint, rayVec);

            % eq 11
            % Calculate the distance from source to target point.
            d12 = vecnorm(rayVec, 2, 2);

            % eq 10
            % Calculate the voxel intersection length.
            tmpDiff = diff(alphas, 1, 2);

            l = d12 .* tmpDiff;
            alphas_mid = alphas(:, 1:end - 1) + 0.5 * tmpDiff;

            % eq 12
            % Calculate the voxel indices: first convert to physical coords
            % and convert to voxel indices
            sourcePoint = matRad_world2cubeCoords(sourcePoint, this.grid, true);
            i = round((sourcePoint(:, 1) + alphas_mid .* rayVec(:, 1)) ./ this.grid.resolution.x);
            j = round((sourcePoint(:, 2) + alphas_mid .* rayVec(:, 2)) ./ this.grid.resolution.y);
            k = round((sourcePoint(:, 3) + alphas_mid .* rayVec(:, 3)) ./ this.grid.resolution.z);

            % Handle numerical instabilities at the borders.
            i(i < 1) = 1;
            i(i > this.numPlanes(1) - 1) = this.numPlanes(1) - 1;
            j(j < 1) = 1;
            j(j > this.numPlanes(2) - 1) = this.numPlanes(2) - 1;
            k(k < 1) = 1;
            k(k > this.numPlanes(3) - 1) = this.numPlanes(3) - 1;

            valIx = ~isnan(alphas_mid);

            % In Matlab direct assignment with sub2ind would work, Octave
            % however does not like NaN values in the subscripts
            ix = NaN(size(valIx));
            ix(valIx) = sub2ind([this.numPlanes(2), this.numPlanes(1), this.numPlanes(3)] - 1, j(valIx), i(valIx), k(valIx));

            for i = 1:numel(this.cubes)
                rho{i} = NaN(size(valIx), class(alphas));
                rho{i}(valIx) = this.cubes{i}(ix(valIx));
            end

        end

    end

    methods (Access = protected)

        function alphas = computeAllAlphas(this, sourcePoint, rayVec)

            % Here we setup grids to enable logical indexing when computing
            % the alphas along each dimension. All alphas between the
            % minimum and maximum index will be computed, with additional
            % exclusion of singular plane occurrences (max == min)
            % All values out of scope will be set to NaN.

            nRays = size(rayVec, 1);

            % eq 4
            % Calculate parametrics values of \alpha_{min} and \alpha_{max} for every
            % axis, intersecting the ray with the sides of the CT.
            aX_1 = (this.planes.x(1) - sourcePoint(:, 1)) ./ rayVec(:, 1);
            aX_end = (this.planes.x(end) - sourcePoint(:, 1)) ./ rayVec(:, 1);

            tmpIx = rayVec(:, 1) == 0;
            aX_1(tmpIx) = NaN;
            aX_end(tmpIx) = NaN;

            aY_1 = (this.planes.y(1) - sourcePoint(:, 2)) ./ rayVec(:, 2);
            aY_end = (this.planes.y(end) - sourcePoint(:, 2)) ./ rayVec(:, 2);

            tmpIx = rayVec(:, 2) == 0;
            aY_1(tmpIx) = NaN;
            aY_end(tmpIx) = NaN;

            aZ_1 = (this.planes.z(1) - sourcePoint(:, 3)) ./ rayVec(:, 3);
            aZ_end = (this.planes.z(end) - sourcePoint(:, 3)) ./ rayVec(:, 3);

            tmpIx = rayVec(:, 3) == 0;
            aZ_1(tmpIx) = NaN;
            aZ_end(tmpIx) = NaN;

            % eq 5
            % Compute the \alpha_{min} and \alpha_{max} in terms of parametric values
            % given by equation 4.
            alphaLimits(:, 1) = max([zeros(nRays, 1) min(aX_1, aX_end) min(aY_1, aY_end) min(aZ_1, aZ_end)], [], 2);
            alphaLimits(:, 2) = min([ones(nRays, 1) max(aX_1, aX_end) max(aY_1, aY_end) max(aZ_1, aZ_end)], [], 2);

            % eq 6
            % Calculate the range of indices who gives parametric values for
            % intersected planes.
            [dimMin, dimMax] = this.computeEntryAndExit(sourcePoint, rayVec, alphaLimits);

            % eq 7
            % For the given range of indices, calculate the paremetrics values who
            % represents intersections of the ray with the plane.
            alphas = this.computePlaneAlphas(sourcePoint, rayVec, dimMin, dimMax);

            % eq 8
            % Merge parametrics sets.
            % The following might look slow but is quite close to Matlab's
            % "unique" implementation
            alphas = sort(horzcat(alphaLimits, alphas{:}), 2); % NaN's are placed at the end when sorting in ascending order
            alphas(diff(alphas, 1, 2) == 0) = NaN; % Remove duplicates
            alphas = sort(alphas, 2); % Again place NaN's at the end

            % Size Reduction (reduce NaN padding) for further computations
            maxNumColumns = max(sum(~isnan(alphas), 2));
            alphas = alphas(:, 1:maxNumColumns);
        end

        function [dimMin, dimMax] = computeEntryAndExit(this, sourcePoint, rayVec, alphaLimits)
            % eq 6
            % Calculate the range of indices who gives parametric values for
            % intersected planes.

            rayDirectionPositive = rayVec > 0;
            alphaLimitsReverse = flip(alphaLimits, 2);

            alphaLimitsRep = repmat(reshape(alphaLimits, [size(alphaLimits, 1) 1 2]), 1, 3, 1);
            alphaLimitsReverseRep = repmat(reshape(alphaLimitsReverse, [size(alphaLimits, 1) 1 2]), 1, 3, 1);
            rayDirPosRep = repmat(rayDirectionPositive, 1, 1, 2);
            alphaAxis = alphaLimitsReverseRep;
            alphaAxis(rayDirPosRep) = alphaLimitsRep(rayDirPosRep);

            tmp = 'xyz';
            [lowerPlanes, upperPlanes] = arrayfun(@(x) deal(this.planes.(x)(1), this.planes.(x)(end)), tmp);
            resArray = arrayfun(@(x) this.grid.resolution.(x), tmp);

            dimMin = this.numPlanes - (upperPlanes - alphaAxis(:, :, 1) .* rayVec - sourcePoint) ./ resArray;
            dimMax = 1 + (sourcePoint + alphaAxis(:, :, 2) .* rayVec - lowerPlanes) ./ resArray;

            dimMin = ceil(round(1e3 * dimMin) / 1e3);
            dimMax = floor(round(1e3 * dimMax) / 1e3);
        end

        function alphas = computePlaneAlphas(this, sourcePoint, rayVec, dimMin, dimMax)
            tmp = 'xyz';
            for i = 1:3
                planeIx = 1:length(this.planes.(tmp(i)));
                planeGrid = repmat(this.planes.(tmp(i)), size(rayVec, 1), 1);
                planeGrid(planeIx < dimMin(:, i) | planeIx > dimMax(:, i) | ...
                          (planeIx == dimMin(:, i) & planeIx == dimMax(:, i)) | ...
                          isnan(dimMin(:, i)) | isnan(dimMax(:, i))) = NaN;
                alphas{i} = (planeGrid - sourcePoint(:, i)) ./ (rayVec(:, i));
            end
        end

    end
end
