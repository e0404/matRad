function coord = matRad_cube2worldCoords(vCoord, ct)
% matRad function to convert cube coordinates to world coordinates
% 
% call
%   coord = matRad_worldToCubeCoordinates(vCoord, ct)
%   
%   ct: matRad ct struct
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
    if isfield(ct, 'dicomInfo') && isfield(ct.dicomInfo,'ImagePositionPatient')
        firstVox = ct.dicomInfo.ImagePositionPatient;
    else 
        firstVox = -(ct.cubeDim./2).*[ct.resolution.x ct.resolution.y ct.resolution.z] ;
    end 
    if prod(vCoord<=ct.cubeDim)  && prod(vCoord>0)       
        coord(1,:) = firstVox(1) + (vCoord(:,1) - 1 ) *ct.resolution.x;
        coord(2,:) = firstVox(2) + (vCoord(:,2) - 1 ) *ct.resolution.y;
        coord(3,:) = firstVox(3) + (vCoord(:,3) - 1 ) *ct.resolution.z;
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