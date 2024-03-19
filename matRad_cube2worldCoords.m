function coord = matRad_cube2worldCoords(vCoord, ct)
% matRad function to convert cube coordinates to world coordinates
% 
% call
%   coord = matRad_worldToCubeCoordinates(vCoord, ct)
%   
% 
%   vCoord : Voxel Coordinates [vx vy vz]
% 
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
if nargin == 2
    if isfield(ct, 'dicomInfo') && isfield(ct.dicomInfo,'ImagePositionPatient')
        firstVox = ct.dicomInfo.ImagePositionPatient;
    else 
        firstVox = - (ct.cubeDim./2).*[ct.resolution.x ct.resolution.y ct.resolution.z] ;
    end 
    try
        coord(1,:) = firstVox(1) + (vCoord(:,1) - 1 ) *ct.resolution.x;
        coord(2,:) = firstVox(2) + (vCoord(:,2) - 1 ) *ct.resolution.y;
        coord(3,:) = firstVox(3) + (vCoord(:,3) - 1 ) *ct.resolution.z;
        coord = coord';
    catch

    end
end

end

