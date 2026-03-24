function [alphas, l, rho, d12, ix] = matRad_siddonRayTracer(isocenterCube, ...
                                                            resolution, ...
                                                            sourcePoint, ...
                                                            targetPoint, ...
                                                            cubes)
% siddon ray tracing through 3D cube to calculate the radiological depth
% according to Siddon 1985 Medical Physics. The raytracer expects the
% isocenter in cube coordinates!
%
% call:
%   [alphas,l,rho,d12,vis] = matRad_siddonRayTracer(isocenter, ...
%                               resolution, ...
%                               sourcePoint, ...
%                               targetPoint, ...
%                               cubes)
%
% input:
%   isocenterCube:  isocenter in cube coordinates [mm]
%   resolution:     resolution of the cubes [mm/voxel]
%   sourcePoint:    source point of ray tracing
%   targetPoint:    target point of ray tracing
%   cubes:          cell array of cubes for ray tracing (it is possible to pass
%                   multiple cubes for ray tracing to save computation time)
%
% output:
%   alphas:         relative distance between start and endpoint for the intersections with the cube
%   l:              lengths of intersections with cubes
%   rho:            densities extracted from cubes
%   d12:            distance between start and endpoint of ray tracing
%   ix:             indices of hit voxels
%
% References
%   [1] http://www.ncbi.nlm.nih.gov/pubmed/4000088
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();
matRad_cfg.dispDeprecationWarning('Calls to matRad_siddonRayTracer will be replaced by usage of the matRad_RayTracerSiddon class in the future!');

grid.resolution = resolution;
grid.dimensions = size(cubes{1});

% At the moment we only implement the Siddon ray-tracer in [1]
hTracer = matRad_RayTracerSiddon(cubes, grid);

isocenterCube = matRad_cubeCoords2worldCoords(isocenterCube, grid);
[alphas, l, rho, d12, ix] = hTracer.traceRay(isocenterCube, sourcePoint, targetPoint);
