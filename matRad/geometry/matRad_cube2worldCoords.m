function coord = matRad_cube2worldCoords(vCoord, gridStruct)
% matRad function to convert cube coordinates to world coordinates
% 
% call
%   coord = matRad_worldToCubeCoordinates(vCoord, gridStruct)
%   
%   gridStruct: matRad ct struct or dij.doseGrid/ctGrid struct
%   vCoord : Voxel Coordinates [vx vy vz]
%   coord: worldCoordinates [x y z ] mm
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coord = [];   % [x y z ] mm 

if nargin == 2 && ~isempty(vCoord)
    if isfield(gridStruct,'cubeDim')
        gridStruct.dimensions = gridStruct.cubeDim;
    end

    if isfield(gridStruct, 'dicomInfo') && isfield(gridStruct.dicomInfo,'ImagePositionPatient')
        firstVox = gridStruct.dicomInfo.ImagePositionPatient;
    else 
        firstVox = -(gridStruct.cubeDim./2).*[gridStruct.resolution.x gridStruct.resolution.y gridStruct.resolution.z] ;
    end 
    if  prod(prod(vCoord<=gridStruct.dimensions))  && prod(prod(vCoord>=0))       
        coord(1,:) = firstVox(1) + (vCoord(:,1) - 1 ) *gridStruct.resolution.x;
        coord(2,:) = firstVox(2) + (vCoord(:,2) - 1 ) *gridStruct.resolution.y;
        coord(3,:) = firstVox(3) + (vCoord(:,3) - 1 ) *gridStruct.resolution.z;
        coord = coord';
    else 
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('Queried cube coordinate is outside the ct cube');
    end
    
else
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Cannot compute coordinates without matRad input coordinates or ct structure');
end

end