function isSame = matRad_isSameGrid(gridA, gridB)
% matRad function to check whether two grids have the same voxel positions
%
% Comparing the dimensions alone is not sufficient: a dose grid whose
% resolution is slightly coarser than the CT resolution can have the same
% number of voxels as the CT, but different voxel positions (e.g. 22 voxels at
% 3 mm and 22 voxels at 3.1 mm, see matRad_createDoseGridAxes). Hence the
% voxel center axes x, y and z are compared as well. If a grid has no axes,
% the resolutions are compared instead.
%
% call:
%   isSame = matRad_isSameGrid(gridA, gridB)
%
% input:
%   gridA, gridB:   grid structs (e.g. dij.ctGrid and dij.doseGrid) with the
%                   field dimensions ([y x z]) and the voxel center axes
%                   x, y, z in mm and/or the resolution.x/y/z in mm
%
% output:
%   isSame:         true if both grids have the same dimensions and voxel
%                   positions (up to 1e-6 mm)
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

tolerance = 1e-6;

isSame = isequal(gridA.dimensions(:), gridB.dimensions(:));

axisNames = {'x', 'y', 'z'};
for i = 1:numel(axisNames)
    if ~isSame
        return
    end
    axisName = axisNames{i};
    if isfield(gridA, axisName) && isfield(gridB, axisName)
        axisA = gridA.(axisName)(:);
        axisB = gridB.(axisName)(:);
        isSame = numel(axisA) == numel(axisB) && all(abs(axisA - axisB) <= tolerance);
    elseif isfield(gridA, 'resolution') && isfield(gridB, 'resolution')
        isSame = abs(gridA.resolution.(axisName) - gridB.resolution.(axisName)) <= tolerance;
    end
end

end
