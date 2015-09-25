function [ix,radDepths,geoDists,x_latDists,z_latDists] = ...
          matRad_calcRadGeoDists(ct, ...
                                 V, ...
                                 isocenter, ...
                                 rot_coords_bev, ...
                                 resolution, ...
                                 sourcePoint, ...
                                 targetPoint, ...
                                 sourcePoint_bev, ...
                                 targetPoint_bev, ...
                                 coords, ...
                                 lateralCutOff, ...
                                 visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad calculation of radiological and geometrical distances used for
% dose calcultion
% 
% call
%   [ix,radDepths,geoDists,x_latDists,z_latDists] = ...
%           matRad_calcRadGeoDists(ct, ...
%                                  V, ...
%                                  isocenter, ...
%                                  rot_coords_bev, ...
%                                  resolution, ...
%                                  sourcePoint, ...
%                                  targetPoint, ...
%                                  sourcePoint_bev, ...
%                                  targetPoint_bev, ...
%                                  coords, ...
%                                  lateralCutOff, ...
%                                  visBool)
%
% input
%   ct:                 ct cube
%   V:                  linear indices of voxels inside the patient, i.e.
%                       potentially interesting voxels
%   isocenter:          isocenter
%   rot_coords_bev:     coordinates of the voxels with index V rotated 
%                       into bev according to the couch and gantry angle        
%   resolution:         resolution of the ct cube [mm]
%   sourcePoint:        source point in voxel coordinates
%   targetPoint:        target point in voxel coordinates
%   sourcePoint_bev:    source point in voxel coordinates in beam's eye view
%   targetPoint_bev:    target point in voxel coordinated in beam's eye view
%   coords:             coordinates of the voxels with index V in standard
%                       coordinate system
%   lateralCutOff:      lateral cutoff specifying the neighbourhood for
%                       which dose calculations will actually be performed
%   visBool:            toogle on/off visualization
%
% output
%   ix:                 indices of voxels where we want to compute dose
%                       influence data
%   radDepths:          corresponding radiological depths
%   geoDists:           corresponding geometrical distances from the
%                       virtual radiation source
%   x_latDists:         lateral x-distance to the central ray (where the
%                       actual computation of the radiological depth takes place)
%   z_latDists:         lateral z-distance to the central ray (where the
%                       actual computation of the radiological depth takes place)
%
% References
%   [1] http://www.ncbi.nlm.nih.gov/pubmed/4000088
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

if nargin < 12
    visBool = 0;
end

% calculate radiological depths on central ray with siddon ray tracer
[alphas,l,rho,d12] = matRad_siddonRayTracer(isocenter,resolution, ...
                                                sourcePoint,targetPoint, ...
                                                {ct});

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

% calculate geometrical distances 
geoDists = ( (coords(ix,1)-sourcePoint(1))*(targetPoint(1) - sourcePoint(1)) + ...
             (coords(ix,2)-sourcePoint(2))*(targetPoint(2) - sourcePoint(2)) + ...
             (coords(ix,3)-sourcePoint(3))*(targetPoint(3) - sourcePoint(3)) ) ...
              / norm(targetPoint-sourcePoint);

% eq 14
% It multiply voxel intersections with \rho values.
% The zero it is neccessary for stability purpose.
d = [0 l .* rho{1}]; %Note. It is not a number "one"; it is the letter "l"

% Calculate accumulated d sum.
dCum = cumsum(d);

% This is necessary for numerical stability.
dCumIx = min([find(dCum==0,1,'last') numel(dCum)-1]);

% Calculate the radiological path
radDepths = interp1(alphas(dCumIx:end),dCum(dCumIx:end),geoDists/d12,'linear',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% visualization
if visBool == 1
        
    x_grid = resolution(1)*[1:size(ct,2)];
    y_grid = resolution(2)*[1:size(ct,1)];
    z_grid = resolution(3)*[1:size(ct,3)];
    [x,y,z] = meshgrid(x_grid,y_grid,z_grid);
    xslice = mean(coords(ix,1)) + isocenter(1);
	yslice = mean(coords(ix,2)) + isocenter(2);
	zslice = mean(coords(ix,3)) + isocenter(3);
    
    figure
    
    subplot(2,3,1)
    colormap bone;
    hold on
    borders = slice(x,y,z,ct,xslice,yslice,zslice);
    set(borders,'FaceColor','interp',...
        'EdgeColor','none',...
        'DiffuseStrength',.8)
    view(-200,27);
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    title('ct')
    
    plot3([sourcePoint(1) targetPoint(1)]+isocenter(1), ...
          [sourcePoint(2) targetPoint(2)]+isocenter(2), ...
          [sourcePoint(3) targetPoint(3)]+isocenter(3),'k')
    
    axis([0 resolution(1)*size(ct,1) ...
          0 resolution(2)*size(ct,2) ...
          0 resolution(3)*size(ct,3)])
    
    colorbar
    box
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % subplot 2
    subplot(2,3,2)
    hold on
    plot(alphas*d12,cumsum(d))
    xlabel('geometric depth along ray [mm]')
    ylabel('radiological depth [mm]')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(2,3,3)
    radDepthCube = zeros(size(ct));
    radDepthCube(V(ix)) = radDepths;
    colormap jet;
    hold on
    borders = slice(x,y,z,radDepthCube,xslice,yslice,zslice);
    set(borders,'FaceColor','interp',...
        'EdgeColor','none',...
        'DiffuseStrength',.8)
    view(-200,27);
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    title('radiological depths')
    
    plot3([sourcePoint(1) targetPoint(1)]+isocenter(1), ...
          [sourcePoint(2) targetPoint(2)]+isocenter(2), ...
          [sourcePoint(3) targetPoint(3)]+isocenter(3),'k')
    
    axis([0 resolution(1)*size(ct,1) ...
          0 resolution(2)*size(ct,2) ...
          0 resolution(3)*size(ct,3)])
      
    colorbar
    box
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(2,3,4)
    geoDistCube = zeros(size(ct));
    geoDistCube(V(ix)) = geoDists;
    colormap jet;
    hold on
    borders = slice(x,y,z,geoDistCube,xslice,yslice,zslice);
    set(borders,'FaceColor','interp',...
        'EdgeColor','none',...
        'DiffuseStrength',.8)
    view(-200,27);
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    title('geometrical distance')
    
    plot3([sourcePoint(1) targetPoint(1)]+isocenter(1), ...
          [sourcePoint(2) targetPoint(2)]+isocenter(2), ...
          [sourcePoint(3) targetPoint(3)]+isocenter(3),'k')
    
    axis([0 resolution(1)*size(ct,1) ...
          0 resolution(2)*size(ct,2) ...
          0 resolution(3)*size(ct,3)])
      
    colorbar
    box
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(2,3,5)
    latDistXCube = zeros(size(ct));
    latDistXCube(V(ix)) = x_latDists;
    colormap jet;
    hold on
    borders = slice(x,y,z,latDistXCube,xslice,yslice,zslice);
    set(borders,'FaceColor','interp',...
        'EdgeColor','none',...
        'DiffuseStrength',.8)
    view(-200,27);
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    title('lateral distance X')
        
    plot3([sourcePoint(1) targetPoint(1)]+isocenter(1), ...
          [sourcePoint(2) targetPoint(2)]+isocenter(2), ...
          [sourcePoint(3) targetPoint(3)]+isocenter(3),'k')
    
    axis([0 resolution(1)*size(ct,1) ...
          0 resolution(2)*size(ct,2) ...
          0 resolution(3)*size(ct,3)])
    
    colorbar
    box
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(2,3,6)
    latDistZCube = zeros(size(ct));
    latDistZCube(V(ix)) = z_latDists;
    colormap jet;
    hold on
    borders = slice(x,y,z,latDistZCube,xslice,yslice,zslice);
    set(borders,'FaceColor','interp',...
        'EdgeColor','none',...
        'DiffuseStrength',.8)
    view(-200,27);
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    title('lateral distance Z')
    
    plot3([sourcePoint(1) targetPoint(1)]+isocenter(1), ...
          [sourcePoint(2) targetPoint(2)]+isocenter(2), ...
          [sourcePoint(3) targetPoint(3)]+isocenter(3),'k')
    
    axis([0 resolution(1)*size(ct,1) ...
          0 resolution(2)*size(ct,2) ...
          0 resolution(3)*size(ct,3)])
      
    colorbar
    box
    
    pause;
    
end