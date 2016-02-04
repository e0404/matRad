function [ix,x_latDists,z_latDists] = ...
          matRad_calcGeoDists(rot_coords_bev, ...
                              sourcePoint_bev, ...
                              targetPoint_bev, ...
                              lateralCutOff)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad calculation of lateral distances from central ray used for
% dose calcultion
% 
% call
%   [ix,x_latDists,z_latDists] = ...
%           matRad_calcRadGeoDists(rot_coords_bev, ...
%                                  sourcePoint_bev, ...
%                                  targetPoint_bev, ...
%                                  lateralCutOff)
%
% input
%   rot_coords_bev:     coordinates of the voxels with index V rotated 
%                       into bev according to the couch and gantry angle        
%   sourcePoint_bev:    source point in voxel coordinates in beam's eye view
%   targetPoint_bev:    target point in voxel coordinated in beam's eye view
%   lateralCutOff:      lateral cutoff specifying the neighbourhood for
%                       which dose calculations will actually be performed
%
% output
%   ix:                 indices of voxels where we want to compute dose
%                       influence data
%   x_latDists:         lateral x-distance to the central ray (where the
%                       actual computation of the radiological depth takes place)
%   z_latDists:         lateral z-distance to the central ray (where the
%                       actual computation of the radiological depth takes place)
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROTATE A SINGLE BEAMLET AND ALIGN WITH BEAMLET WHO PASSES THROUGH
% ISOCENTER

% Put [0 0 0] position in the source point for beamlet who passes through
% isocenter
a = (-sourcePoint_bev - sourcePoint_bev)';

% Normalize the vector
a = a/norm(a);

% Put [0 0 0] position in the source point for a single beamlet
b = (targetPoint_bev - sourcePoint_bev)';

% Normalize the vector
b = b/norm(b);

% Define function for obtain rotation matrix.
if sum(a==b)==3 % rotation matrix corresponds to eye matrix if the vectors are the same
    RU = @(a,b) eye(3);
else
    % Define fuction to obtain skew symmetric cross-product matrix of vector v
    ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
    RU = @(a,b) eye(3) + ssc(cross(a,b)) + ssc(cross(a,b))^2*(1-dot(a,b))/(norm(cross(a,b))^2);
end

% Calculate rotation matrix for rotate a single beamlet to be aligned to
% beamlet who passes through isocenter.
R = RU(a,b);

% Rotate every CT voxel 
rot_coords_temp = rot_coords_bev*R;

% Put [0 0 0] position CT in center of the beamlet.
x_latDists = rot_coords_temp(:,1) + sourcePoint_bev(1);
z_latDists = rot_coords_temp(:,3) + sourcePoint_bev(3);

rad_distancesSq = x_latDists.^2 + z_latDists.^2;
ix = find(rad_distancesSq <= lateralCutOff^2);

x_latDists = x_latDists(ix);
z_latDists = z_latDists(ix);
