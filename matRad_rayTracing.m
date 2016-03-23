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

% set up rad depth cube for results
radDepthCube = NaN*ones(size(ct.cube));

% set up coordinates of all voxels in cube
[xCoords_vox, yCoords_vox, zCoords_vox] = meshgrid(1:size(ct.cube,1),1:size(ct.cube,2),1:size(ct.cube,3));

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
rayMx_bev_y = max(coords_bev(V,2)) + max([ct.resolution.x ct.resolution.y ct.resolution.z]);


xCoords = xCoords-stf.sourcePoint(1);
yCoords = yCoords-stf.sourcePoint(2);
zCoords = zCoords-stf.sourcePoint(3);
coords  = [xCoords yCoords zCoords];
    
% calculate geometric distances
if nargout > 1
    geoDistCube = reshape(sqrt(sum(coords.^2,2)),size(ct.cube));
end

% calculate spacing of rays on ray matrix
rayMxSpacing = min([ct.resolution.x ct.resolution.y ct.resolution.z]);

% define candidate ray matrix covering 1000x1000mm^2
numOfCandidateRays = 2 * ceil(500/rayMxSpacing) + 1;
candidateRayMx     = zeros(numOfCandidateRays);

% define coordinates
[candidateRaysCoords_X,candidateRaysCoords_Z] = meshgrid(rayMxSpacing*[floor(-500/rayMxSpacing):ceil(500/rayMxSpacing)]);

% check which rays should be used
for i = 1:stf.numOfRays
   
    ix = (candidateRaysCoords_X(:)-stf.ray(i).rayPos_bev(1)).^2 + ...
         (candidateRaysCoords_Z(:)-stf.ray(i).rayPos_bev(3)).^2 ...
           <= lateralCutoff^2;
    
    candidateRayMx(ix) = 1;
    
end

% set up ray matrix
rayMx_bev = [candidateRaysCoords_X(logical(candidateRayMx(:))) ...
              rayMx_bev_y*ones(sum(candidateRayMx(:)),1) ...  
              candidateRaysCoords_Z(logical(candidateRayMx(:)))];

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
metricHitVoxelsCube = -inf*ones(size(ct.cube));

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

    dotProdHitVoxels = coords(ixHitVoxel,:)*normRayVector';
    
    ixRememberFromCurrTracing = dotProdHitVoxels > metricHitVoxelsCube(ixHitVoxel)';

    if any(ixRememberFromCurrTracing) > 0
        
        metricHitVoxelsCube(ixHitVoxel(ixRememberFromCurrTracing)) = dotProdHitVoxels(ixRememberFromCurrTracing);

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


