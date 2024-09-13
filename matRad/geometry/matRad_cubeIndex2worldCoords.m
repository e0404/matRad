function coord = matRad_cubeIndex2worldCoords(cubeIx, gridStruct)
% matRad function to convert cube indices to world coordinates
% 
% call
%   coord = matRad_cubeIndex2worldCoords(vCoord, gridStruct)
% 
% inputs
%   cCoord:         cube indices [i j k] (Nx3) or [linIx] (Nx1)
%   gridStruct:     matRad ct struct or dij.doseGrid/ctGrid struct
%   allowOutside:   indices not within the image bounds will be calculated
%                   optional, default is true
%
% outputs
%   coord:          worldCoordinates [x y z] (Nx3 in mm)
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sanitize grid / cube dimensions
if isfield(gridStruct,'cubeDim')
    gridStruct.dimensions = gridStruct.cubeDim;
end

%Check if we have linear indices
if size(cubeIx,2) == 1
    [cubeIx(:,1),cubeIx(:,2),cubeIx(:,3)] = ind2sub(gridStruct.dimensions,cubeIx);
end

%Check if we have the right dimensions
if size(cubeIx,2) ~= 3
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('voxel coordinates must be Nx3 (subscript indices) or Nx1 (linear indices)!');
end

%First create cube coordinates
coord = cubeIx(:,[2 1 3]) .* [gridStruct.resolution.x gridStruct.resolution.y gridStruct.resolution.z];
coord = matRad_cubeCoords2worldCoords(coord,gridStruct,false);
    
end