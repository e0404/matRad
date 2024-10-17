function coord = matRad_world2cubeCoords(wCoord, gridStruct , allowOutside)
% matRad function to convert world coordinates to cube coordinates
% 
% call
%   coord = world2cubeCoords(wCoord, ct)
%
%
%   wCoord:         world coordinates array Nx3 (x,y,z) [mm]
%   gridStruct:     can be matRad ct, dij.doseGrid, or the ctGrid
%                   required fields x,y,x,dimensions,resolution
%   allowOutside:   indices not within the image bounds will be calculated
%                   optional, default is true
%
%   coord :         cube coordinates (x,y,z) - Nx3 in mm    
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
translation   = firstVoxCube - firstVoxWorld;

% If we don't allow outside coordinates, we check here
if (~allowOutside)
    % world bounds
    minBoundsWorld = [min(gridStruct.x)-gridStruct.resolution.x/2 min(gridStruct.y)-gridStruct.resolution.y/2 min(gridStruct.z)-gridStruct.resolution.z/2];
    maxBoundsWorld = [max(gridStruct.x)+gridStruct.resolution.x/2 max(gridStruct.y)+gridStruct.resolution.y/2 max(gridStruct.z)+gridStruct.resolution.z/2];
    violateBounds = any(any(wCoord < minBoundsWorld | wCoord > maxBoundsWorld));

    if violateBounds
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('Queried world coordinate is outside the ct');
    end
end
    
% find the closest world coord index in gridStruct.x/y/z
% calc absolute differences and locate smallest difference
coord = wCoord + translation;

end
