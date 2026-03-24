function [radDepthV, radDepthCube] = matRad_rayTracing(stfElement, ct, V, rot_coordsV, lateralCutoff)
% matRad visualization of two-dimensional dose distributions on ct including
% segmentation
%
% call
%   [radDepthV, radDepthCube] = matRad_rayTracing(stf,ct,V,rot_coordsV,lateralCutoff)
%
% input
%   stfElement:    matRad steering information struct of single(!) beam
%   ct:            ct cube
%   V:             linear voxel indices e.g. of voxels inside patient.
%   rot_coordsV    coordinates in beams eye view inside the patient
%   lateralCutoff: lateral cut off used for ray tracing (optional)

%
% output
%   radDepthV:      radiological depth inside the patient
%   radDepthCube:   radiological depth in whole ct
%
% References
%   [1] http://www.sciencedirect.com/science/article/pii/S1120179711001359
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
matRad_cfg.dispDeprecationWarning('Calls to matRad_rayTracing will be replaced by usage of the matRad_RayTracer class in the future!');

% At the moment we only implement the Siddon ray-tracer in [1]
hTracer = matRad_RayTracerSiddon(ct.cube, ct);

if nargin >= 5
    hTracer.lateralCutOff = lateralCutoff;
end

if nargin >= 4
    [radDepthV, radDepthCube] = hTracer.traceCube(stfElement, V, rot_coordsV);
elseif nargin >= 3
    [radDepthV, radDepthCube] = hTracer.traceCube(stfElement, V);
else
    [radDepthV, radDepthCube] = hTracer.traceCube(stfElement);
end
