function [ix,rad_distancesSq,isoLatDistsX,isoLatDistsZ,latDistsX,latDistsZ] = ...
          matRad_calcGeoDistsXXXfast(rot_coords_bev, ...
                              sourcePoint_bev, ...
                              targetPoint_bev, ...
                              SAD, ...
                              radDepthIx, ...
                              lateralCutOff)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad calculation of lateral distances from central ray used for
% dose calcultion
% 
% call
%   [ix,x_latDists,z_latDists] = ...
%           matRad_calcGeoDists(rot_coords_bev, ...
%                               sourcePoint_bev, ...
%                               targetPoint_bev, ...
%                               SAD, ...
%                               radDepthIx, ...
%                               lateralCutOff)
%
% input
%   rot_coords_bev:     coordinates in bev of the voxels with index V,
%                       where also ray tracing results are availabe 
%   sourcePoint_bev:    source point in voxel coordinates in beam's eye view
%   targetPoint_bev:    target point in voxel coordinated in beam's eye view
%   SAD:                source-to-axis distance
%   radDepthIx:         sub set of voxels for which radiological depth
%                       calculations are available
%   lateralCutOff:      lateral cutoff specifying the neighbourhood for
%                       which dose calculations will actually be performed
%
% output
%   ix:                 indices of voxels where we want to compute dose
%                       influence data
%   isoLatDistsX:       lateral x-distance to the central ray projected to
%                       iso center plane
%   isoLatDistsZ:       lateral z-distance to the central ray projected to
%                       iso center plane
%   radialDist_sq:      squared radial distance to the central ray (where the
%                       actual computation of the radiological depth takes place)
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROTATE A SINGLE BEAMLET AND ALIGN WITH BEAMLET WHO PASSES THROUGH
% ISOCENTER

% Define function for obtain rotation matrix.
rot_coords_temp = rot_coords_bev;

% Put [0 0 0] position CT in center of the beamlet.
latDistsX = rot_coords_temp(:,1) + sourcePoint_bev(1);
latDistsZ = rot_coords_temp(:,3) + sourcePoint_bev(3);

% check of radial distance exceeds lateral cutoff (projected to iso center)
rad_distancesSq = latDistsX.^2 + latDistsZ.^2;
subsetMask = rad_distancesSq ./ rot_coords_temp(:,2).^2 <= lateralCutOff^2 /SAD^2;

% return index list within considered voxels
ix = radDepthIx(subsetMask);

% return radial distances squared
if nargout > 1
    rad_distancesSq = rad_distancesSq(subsetMask);
end

% return x & z distance
if nargout > 2
   isoLatDistsX = latDistsX(subsetMask)./rot_coords_temp(subsetMask,2)*SAD;
   isoLatDistsZ = latDistsZ(subsetMask)./rot_coords_temp(subsetMask,2)*SAD; 
end

latDistsX = latDistsX(subsetMask);
latDistsZ = latDistsZ(subsetMask);