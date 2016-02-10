function [radDepthCube,geoDistCube] = matRad_rayTracing(stf,ct,V,lateralCutoff)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad visualization of two-dimensional dose distributions on ct including
% segmentation
% 
% call
%   [radDepthCube,geoDistCube] = matRad_rayTracing(stf,ct,V,lateralCutoff)
%
% input
%   stf:           matRad steering information struct of one beam
%   ct:            ct cube
%   V:             linear voxel indices e.g. of voxels inside patient.
%   lateralCutoff: lateral cut off used for ray tracing

%
% output
%   radDepthCube:  radiological depth cube in the ct.cube dimensions
%   geoDistCube:   geometrical distance cube in the ct.cube dimensions
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

% set up rad depth cube for results
radDepthCube = inf*ones(size(ct.cube));

% set up coordinates of all voxels in cube
[yCoords_vox, xCoords_vox, zCoords_vox] = ind2sub(size(ct.cube),1:numel(ct.cube));

xCoords = xCoords_vox(:)*ct.resolution.x-stf.isoCenter(1);
yCoords = yCoords_vox(:)*ct.resolution.y-stf.isoCenter(2);
zCoords = zCoords_vox(:)*ct.resolution.z-stf.isoCenter(3);

% Rotation around Z axis (gantry)
inv_rotMx_XY_T = [ cosd(-stf.gantryAngle) sind(-stf.gantryAngle) 0;
                  -sind(-stf.gantryAngle) cosd(-stf.gantryAngle) 0;
                                        0                      0 1];

% Rotation around Y axis (Couch movement)
inv_rotMx_XZ_T = [cosd(-stf.couchAngle) 0 -sind(-stf.couchAngle);
                                      0 1                      0;
                  sind(-stf.couchAngle) 0  cosd(-stf.couchAngle)];
 

coords_bev = [xCoords yCoords zCoords]*inv_rotMx_XZ_T*inv_rotMx_XY_T;             
              
% set up ray matrix direct behind last voxel
rayMx_bev_z = max(coords_bev(V,2)) + 1;


xCoords = xCoords-stf.sourcePoint(1);
yCoords = yCoords-stf.sourcePoint(2);
zCoords = zCoords-stf.sourcePoint(3);
coords  = [xCoords yCoords zCoords];
    
% calculate geometric distances
if nargout > 1
    geoDistCube = sqrt(sum(coords.^2,2));
end
% set up ray matrix
rayMxSpacing = min([ct.resolution.x ct.resolution.y ct.resolution.z]);

rayMx_bev    = [];
numOfRayTracingRays    = ceil((stf.sourcePoint_bev(2)-rayMx_bev_z)/stf.sourcePoint_bev(2) * lateralCutoff / rayMxSpacing);

for j = 1:stf.numOfRays

    tmp_rayMx_bev_x = repmat(round(stf.ray(j).rayPos_bev(1) / rayMxSpacing) + [-numOfRayTracingRays:numOfRayTracingRays],2*numOfRayTracingRays+1,1);
    tmp_rayMx_bev_y = rayMx_bev_z * ones(2*numOfRayTracingRays+1);
    tmp_rayMx_bev_z = repmat(round(stf.ray(j).rayPos_bev(3) / rayMxSpacing) + [-numOfRayTracingRays:numOfRayTracingRays]',1,2*numOfRayTracingRays+1);

    rayMx_bev = [rayMx_bev; [rayMxSpacing * tmp_rayMx_bev_x(:) tmp_rayMx_bev_y(:) rayMxSpacing * tmp_rayMx_bev_z(:)]];

    rayMx_bev = unique(rayMx_bev,'rows');

end

%     figure,
%     for jj = 1:length(rayMx_bev)
%        plot(rayMx_bev(jj,1),rayMx_bev(jj,3),'rx'),hold on 
%     end
    
% Rotation around Z axis (gantry)
rotMx_XY_T = [ cosd(stf.gantryAngle) sind(stf.gantryAngle) 0;
              -sind(stf.gantryAngle) cosd(stf.gantryAngle) 0;
                                   0                     0 1];
    
% Rotation around Y axis (couch)
rotMx_XZ_T = [cosd(stf.couchAngle) 0 -sind(stf.couchAngle);
                                 0 1                     0;
              sind(stf.couchAngle) 0  cosd(stf.couchAngle)];

% rotate ray matrix from bev to world coordinates
rayMx_world = rayMx_bev * rotMx_XY_T * rotMx_XZ_T;

% set up distance cube to decide which rad depths should be stored
rayTracingDotProdCube = -inf*ones(size(ct.cube));

% perform ray tracing over all rays
for j = 1:size(rayMx_world,1)

    % run siddon ray tracing algorithm
    [alphas,l,rho,d12,ixHitVoxel] = matRad_siddonRayTracer(stf.isoCenter, ...
                                ct.resolution, ...
                                stf.sourcePoint, ...
                                rayMx_world(j,:), ...
                                {ct.cube});
                                                        
    % find voxels for which we should remember this tracing because this is
    % the closest ray
    normRayVector = rayMx_world(j,:) - stf.sourcePoint;
    normRayVector = normRayVector/norm(normRayVector);

    dotProdHitVoxels = coords(ixHitVoxel,:)*normRayVector'; % this also corresponds to the geometrical distance!!!

    ixRememberFromCurrTracing = dotProdHitVoxels > rayTracingDotProdCube(ixHitVoxel)';

    if sum(ixRememberFromCurrTracing) > 0
        rayTracingDotProdCube(ixHitVoxel(ixRememberFromCurrTracing)) = dotProdHitVoxels(ixRememberFromCurrTracing);

        % calc radioloical depths

        % eq 14
        % It multiply voxel intersections with \rho values.
        % The zero it is neccessary for stability purpose.
        d = [0 l .* rho{1}]; %Note. It is not a number "one"; it is the letter "l"

        % Calculate accumulated d sum.
        dCum = cumsum(d);

        % Calculate the radiological path
        vRadDepth = interp1(alphas,dCum,dotProdHitVoxels(ixRememberFromCurrTracing)/d12,'linear',0);
             
        radDepthCube(ixHitVoxel(ixRememberFromCurrTracing)) = vRadDepth;
    end
    
end

