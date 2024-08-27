function coord = matRad_cubeCoords2worldCoords(cCoord, gridStruct, allowOutside)
% matRad function to convert cube coordinates to world coordinates
% 
% call
%   coord = matRad_worldToCubeCoordinates(vCoord, gridStruct)
% 
% inputs
%   cCoord:         cube coordinates [vx vy vz] (Nx3 in mm)
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

if nargin < 3
    allowOutside = true;
end
          
gridStruct = matRad_getWorldAxes(gridStruct);  

firstVoxWorld = [min(gridStruct.x)  min(gridStruct.y) min(gridStruct.z)];
firstVoxCube  = [gridStruct.resolution.x gridStruct.resolution.y gridStruct.resolution.z];
translation   = firstVoxWorld - firstVoxCube;

% If we don't allow outside coordinates, we check here
if (~allowOutside)
    minBoundsCube = firstVoxCube./2;
    maxBoundsCube = firstVoxCube./2 + gridStruct.dimensions([2 1 3]).*[gridStruct.resolution.x gridStruct.resolution.y gridStruct.resolution.z];
    violateBounds = any(any(cCoord < minBoundsCube | cCoord > maxBoundsCube));

    if violateBounds
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('Queried world coordinate is outside the ct');
    end
end
    
% find the closest world coord index in gridStruct.x/y/z
% calc absolute differences and locate smallest difference
coord = cCoord + translation;