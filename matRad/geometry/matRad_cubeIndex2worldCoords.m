function coord = matRad_cubeIndex2worldCoords(cubeIx, gridStruct)
% matRad function to convert cube indices to world coordinates
%
% call:
%   coord = matRad_cubeIndex2worldCoords(vCoord, gridStruct)
%
% input:
%   cCoord:         cube indices [i j k] (Nx3) or [linIx] (Nx1)
%   gridStruct:     matRad ct struct or dij.doseGrid/ctGrid struct
%   allowOutside:   indices not within the image bounds will be calculated
%                   optional, default is true
%
% output:
%   coord:          worldCoordinates [x y z] (Nx3 in mm)
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sanitize grid / cube dimensions
if isfield(gridStruct, 'cubeDim')
    gridStruct.dimensions = gridStruct.cubeDim;
end

% Check if we have linear indices
if size(cubeIx, 2) == 1
    [s1, s2, s3] = ind2sub(gridStruct.dimensions, cubeIx);
    cubeIx = [s2, s1, s3]; % ijk/xyz -> jik
elseif size(cubeIx, 2) == 3 % if we have subscript coordinates, permute ijk/xyz -> jik
    cubeIx = cubeIx(:, [2 1 3]);
else % we have the wrong dimensions
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('voxel coordinates must be Nx3 (subscript indices) or Nx1 (linear indices)!');
end

% First create cube coordinates
coord = double(cubeIx) .* [gridStruct.resolution.x gridStruct.resolution.y gridStruct.resolution.z];
coord = matRad_cubeCoords2worldCoords(coord, gridStruct, false);

end
