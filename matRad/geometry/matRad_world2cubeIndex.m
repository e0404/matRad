function indices = matRad_world2cubeIndex(wCoord, gridStruct,allowOutside)
% matRad function to convert world coordinates to cube indices
% 
% call
%   coord = world2cubeCoords(wCoord, ct)
%
%
%   wCoord:         world coordinates array Nx3 (x,y,z) [mm]
%   gridStruct:     can be matRad ct, dij.doseGrid, or the ctGrid
%                   required fields x,y,x,dimensions,resolution
%   allowOutside:   If coordinates outside are allowed. False default.
%
%   index:          cube index (i,j,k) honoring Matlab permuatation of i,j
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
    allowOutside = false;
end

if size(wCoord,2) ~= 3
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('world coordinates must be Nx3');
end    
          
gridStruct = matRad_getWorldAxes(gridStruct);  

%First obtain coordinates in multiples of the resolution
coords = matRad_world2cubeCoords(wCoord,gridStruct,allowOutside);

%Now get the indices
coords = round(coords ./ [gridStruct.resolution.x gridStruct.resolution.y gridStruct.resolution.z]);

%Do the permutation
indices = coords(:,[2 1 3]);

end
