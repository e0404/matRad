function [x, y, z] = matRad_createDoseGridAxes(ctGrid, doseResolution)
% matRad function to create the voxel center axes of a dose grid for a given CT grid
%
% The dose grid is centered on the CT: the midpoint of the dose grid voxel
% centers coincides with the midpoint of the CT voxel centers, and the number
% of dose voxels is the smallest one for which the dose grid covers the
% complete CT volume (the outer voxel faces of the dose grid lie at or beyond
% the outer voxel faces of the CT). Consequently, the outermost dose voxel
% centers may lie up to half a dose voxel outside of the range of the CT voxel
% centers, but never outside of the CT volume. Whenever data is resampled
% between the grids, it needs to be extrapolated accordingly (see
% matRad_interp3 with the extrapolation mode 'nearest').
%
% If the dose resolution along an axis equals the CT resolution, the dose grid
% axis is an exact copy of the CT axis.
%
% call:
%   [x, y, z] = matRad_createDoseGridAxes(ctGrid, doseResolution)
%
% input:
%   ctGrid:         CT grid struct with the world axes x, y, z (voxel
%                   centers in mm) and resolution.x/y/z in mm
%   doseResolution: struct with the dose grid resolution fields x, y, z in mm
%
% output:
%   x, y, z:        row vectors with the voxel center coordinates of the
%                   dose grid in mm
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

x = matRad_createCenteredAxis(ctGrid.x, ctGrid.resolution.x, doseResolution.x);
y = matRad_createCenteredAxis(ctGrid.y, ctGrid.resolution.y, doseResolution.y);
z = matRad_createCenteredAxis(ctGrid.z, ctGrid.resolution.z, doseResolution.z);

end

function doseAxis = matRad_createCenteredAxis(ctAxis, ctRes, doseRes)
% Centered, covering dose axis along a single dimension

% Same resolution: reproduce the CT axis exactly (no rounding of the formula
% below) so that equality checks between the grids hold
if abs(doseRes - ctRes) <= 1e-9 * abs(ctRes)
    doseAxis = ctAxis(:)';
else
    % Tolerance against ratios that are integer mathematically but slightly
    % above in floating point (e.g. 30 * 0.1 / 0.3)
    tolerance = 1e-6;
    numOfVoxels = max(1, ceil(numel(ctAxis) * ctRes / doseRes - tolerance));

    ctCenter = (ctAxis(1) + ctAxis(end)) / 2;
    doseAxis = ctCenter + ((0:numOfVoxels - 1) - (numOfVoxels - 1) / 2) * doseRes;
end

end
