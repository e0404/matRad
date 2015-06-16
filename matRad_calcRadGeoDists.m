function [ix,radDepths,geoDists,x_latDists,z_latDists] = ...
          matRad_calcRadGeoDists(ct,V,isocenter,rot_coords,...
                                 resolution,sourcePoint,targetPoint,sourcePoint_bev,...
                                 targetPoint_bev,X,Y,Z,lateralCutOff,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad calculation of radiological and geometrical distances used for
% dose calcultion
% 
% call
%   [ix,radDepths,geoDists,x_latDists,z_latDists] = ...
%          matRad_calcRadGeoDists(ct,V,isocenter,rot_coords,...
%                                 resolution,sourcePoint,targetPoint,sourcePoint_bev,...
%                                 targetPoint_bev,X,Y,Z,lateralCutOff,visBool)
%
% input
%   ct:                 ct cube
%   V:                  linear indices of voxels inside the patient, i.e.
%                       potentially interesting voxels
%   isocenter:          isocenter
%   rot_coords:         coordinates of the voxels with index V rotated 
%                       according to the couch and gantry angle        
%   resolution:         resolution of the ct cube [mm]
%   sourcePoint:        source point in voxel coordinates
%   targetPoint:        target point in voxel coordinates
%   sourcePoint_bev:    source point in voxel coordinates in beam's eye view
%   targetPoint_bev:    target point in voxel coordinated in beam's eye view
%   X:                  x coordinates of the voxels with index V
%   Y:                  y coordinates of the voxels with index V
%   Z:                  z coordinates of the voxels with index V
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

% calculate radiological depths on central ray with siddon ray tracer
[alphas,l,rho,d12,vis] = matRad_siddonRayTracer(isocenter,resolution, ...
                                                sourcePoint,targetPoint, ...
                                                {ct},visBool);
                                           
% Add isocenter to source and target point. Because the algorithm does not
% works with negatives values. This puts (0,0,0) in the center of first
% voxel
sourcePoint = sourcePoint + isocenter;
targetPoint = targetPoint + isocenter;

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

% Put [0 0 0] position in the source point for every CT voxel. This is
% necessary because rotation is against source point.
rot_coords_bev(:,1) = rot_coords(:,1)-sourcePoint_bev(1);
rot_coords_bev(:,2) = rot_coords(:,2)-sourcePoint_bev(2);
rot_coords_bev(:,3) = rot_coords(:,3)-sourcePoint_bev(3);

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
geoDists = ( (X(ix)-sourcePoint(1))*(targetPoint(1) - sourcePoint(1)) + ...
             (Y(ix)-sourcePoint(2))*(targetPoint(2) - sourcePoint(2)) + ...
             (Z(ix)-sourcePoint(3))*(targetPoint(3) - sourcePoint(3)) ) ...
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
radDepths = interp1(alphas(dCumIx:end)*d12,dCum(dCumIx:end),geoDists,'linear',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% visualization
if visBool == 1
    
    %Save the numbers of planes.
    [xNumPlanes, yNumPlanes, zNumPlanes] = size(ct);
    y_latDists = rot_coords_temp(:,2) + sourcePoint_bev(2);
    y_latDists = y_latDists(ix);
    
    %Data with less than a certain of radial distance from the beamlet.
    ix_vis = rad_distancesSq > lateralCutOff^2;
    V(ix_vis) = [];

    
    close all;
    figure;
    subplot(1, 2, 1);
	plot3(rot_coords(:,1),rot_coords(:,2),rot_coords(:,3),'kx');
    hold on;
    plot3([sourcePoint_bev(1) targetPoint_bev(1)],[sourcePoint_bev(2) targetPoint_bev(2)],[sourcePoint_bev(3) targetPoint_bev(3)],'g');
    plot3([0 0],[1000 -1000],[0 0],'r');
    xlim([-400 400])
    ylim([-1000 1000])
    zlim([-400 400])
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    title('Beamlet unaligned')
    box
    
    subplot(1, 2, 2);
    new_target_bev = (targetPoint_bev-sourcePoint_bev)*R+sourcePoint_bev;
    plot3(x_latDists,y_latDists,z_latDists,'kx');
    hold on;
    plot3([new_target_bev(1) sourcePoint_bev(1)],[new_target_bev(2) sourcePoint_bev(2)],[new_target_bev(3) sourcePoint_bev(3)],'g');
    plot3([0 0],[1000 -1000],[0 0],'r');
    xlim([-400 400])
    ylim([-1000 1000])
    zlim([-400 400])
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    title('Beamlet aligned')
    box
    
    if ~isempty(vis.alpha_x)
        coor_x(:,1) = sourcePoint(1) + vis.alpha_x * (targetPoint(1) - sourcePoint(1));
        coor_x(:,2) = sourcePoint(2) + vis.alpha_x * (targetPoint(2) - sourcePoint(2));
        coor_x(:,3) = sourcePoint(3) + vis.alpha_x * (targetPoint(3) - sourcePoint(3));
        inside_ct_x = coor_x(:,1) < 0 | coor_x(:,1) > xNumPlanes*resolution(1) | coor_x(:,2) < 0 | coor_x(:,2) > yNumPlanes*resolution(2) | coor_x(:,3) < 0 | coor_x(:,3) > zNumPlanes*resolution(3);
        coor_x(inside_ct_x,:) = [];
    else
        coor_x = [];
    end
    
    if ~isempty(vis.alpha_y)
        coor_y(:,1) = sourcePoint(1) + vis.alpha_y * (targetPoint(1) - sourcePoint(1));
        coor_y(:,2) = sourcePoint(2) + vis.alpha_y * (targetPoint(2) - sourcePoint(2));
        coor_y(:,3) = sourcePoint(3) + vis.alpha_y * (targetPoint(3) - sourcePoint(3));
        inside_ct_y = coor_y(:,1) < 0 | coor_y(:,1) > xNumPlanes*resolution(1) | coor_y(:,2) < 0 | coor_y(:,2) > yNumPlanes*resolution(2) | coor_y(:,3) < 0 | coor_y(:,3) > zNumPlanes*resolution(3);
        coor_y(inside_ct_y,:) = [];
    else
        coor_y = [];
    end
    
    if ~isempty(vis.alpha_z)
        coor_z(:,1) = sourcePoint(1) + vis.alpha_z * (targetPoint(1) - sourcePoint(1));
        coor_z(:,2) = sourcePoint(2) + vis.alpha_z * (targetPoint(2) - sourcePoint(2));
        coor_z(:,3) = sourcePoint(3) + vis.alpha_z * (targetPoint(3) - sourcePoint(3));
        inside_ct_z = coor_z(:,1) < 0 | coor_z(:,1) > xNumPlanes*resolution(1) | coor_z(:,2) < 0 | coor_z(:,2) > yNumPlanes*resolution(2) | coor_z(:,3) < 0 | coor_z(:,3) > zNumPlanes*resolution(3);
        coor_z(inside_ct_z,:) = [];
    else
        coor_z = [];
    end
    figure;
    colormap jet;
    clf
    x_min = vis.x_min -isocenter(1);
    y_min = vis.y_min -isocenter(2);
    z_min = vis.z_min -isocenter(3);
    x_max = vis.x_max -isocenter(1);
    y_max = vis.y_max -isocenter(2);
    z_max = vis.z_max -isocenter(3);
    
    subplot(2,3,1)
    colormap bone;
    hold on
    if ~isempty(coor_x)
        plot3(coor_x(:,1)-isocenter(1),coor_x(:,2)-isocenter(2),coor_x(:,3)-isocenter(3),'mx');
    end
    if ~isempty(coor_y)
        plot3(coor_y(:,1)-isocenter(1),coor_y(:,2)-isocenter(2),coor_y(:,3)-isocenter(3),'cx');
    end
    if ~isempty(coor_z)
        plot3(coor_z(:,1)-isocenter(1),coor_z(:,2)-isocenter(2),coor_z(:,3)-isocenter(3),'yx');
    end
    plot3(x_min,y_min,z_min,'or');
    plot3(x_max,y_max,z_max,'og');
    xlim([0*resolution(1)-isocenter(1),xNumPlanes*resolution(1)-isocenter(1)])
    ylim([0*resolution(2)-isocenter(2),yNumPlanes*resolution(2)-isocenter(2)])
    zlim([0*resolution(3)-isocenter(3),zNumPlanes*resolution(3)-isocenter(3)])
    view(-200,45);
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    box
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % subplot 2
    subplot(2,3,2)
    hold on
    plot(alphas*d12,cumsum(d))
    xlabel('geometric depth along ray [mm]')
    ylabel('radiological depth [mm]')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_grid = 0.5*resolution(1)-isocenter(1):resolution(1):xNumPlanes*resolution(1)-isocenter(1);
    y_grid = 0.5*resolution(2)-isocenter(2):resolution(2):yNumPlanes*resolution(2)-isocenter(2);
    z_grid = 0.5*resolution(3)-isocenter(3):resolution(3):zNumPlanes*resolution(3)-isocenter(3);
    xslice = 0;
	yslice = 0;
	zslice = 0;
    
    subplot(2,3,3)
     colormap jet;
    [x,y,z] = meshgrid(x_grid,y_grid,z_grid);
    %radDepths_vis = zeros(size(ct));
    radDepths_vis = zeros(size(ct));
    radDepths_vis(V) = radDepths;
    slice(x,y,z,radDepths_vis,1,1,1) % Draw some volume boundaries
    borders = slice(x,y,z,radDepths_vis,xslice,yslice,zslice);
    set(borders,'FaceColor','interp',...
        'EdgeColor','none',...
        'DiffuseStrength',.8)
    xlim([0.5*resolution(1)-isocenter(1),xNumPlanes*resolution(1)-isocenter(1)])
    ylim([0.5*resolution(2)-isocenter(2),yNumPlanes*resolution(2)-isocenter(2)])
    zlim([0.5*resolution(3)-isocenter(3),zNumPlanes*resolution(3)-isocenter(3)])
    hold on;
    plot3(x_min,y_min,z_min,'or');
    plot3(x_max,y_max,z_max,'og');
    plot3([x_min x_max],[y_min y_max],[z_min z_max],'w')
    hold off
    grid off
    view(-200,27);
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    title('Radiological depths')
    colorbar
    box

    %%%%%%%%%%%%%%
    subplot(2,3,4)
    colormap jet;
    x_latDists_vis = zeros(size(ct));
	x_latDists_vis(V) = x_latDists;
    slice(x,y,z,x_latDists_vis,1,1,1) % Draw some volume boundaries
    borders = slice(x,y,z,x_latDists_vis,xslice,yslice,zslice);
    set(borders,'FaceColor','interp',...
        'EdgeColor','none',...
        'DiffuseStrength',.8)
    xlim([0*resolution(1)-isocenter(1),xNumPlanes*resolution(1)-isocenter(1)])
    ylim([0*resolution(2)-isocenter(2),yNumPlanes*resolution(2)-isocenter(2)])
    zlim([0*resolution(3)-isocenter(3),zNumPlanes*resolution(3)-isocenter(3)])
    hold on;
    plot3(x_min,y_min,z_min,'or');
    plot3(x_max,y_max,z_max,'og');
    plot3([x_min x_max],[y_min y_max],[z_min z_max],'w');
    hold off
    grid off
    view(-200,27);
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    title('X lateral distances')
    colorbar
    box

    %%%%%%%%%%%%%%
    subplot(2,3,5)
	y_latDists_inside = zeros(size(ct));
	y_latDists_inside(V) = y_latDists;
	y_latDists = y_latDists_inside;
    slice(x,y,z,y_latDists,1,1,1) % Draw some volume boundaries
    borders = slice(x,y,z,y_latDists,xslice,yslice,zslice);
    set(borders,'FaceColor','interp',...
        'EdgeColor','none',...
        'DiffuseStrength',.8)
    xlim([0*resolution(1)-isocenter(1),xNumPlanes*resolution(1)-isocenter(1)])
    ylim([0*resolution(2)-isocenter(2),yNumPlanes*resolution(2)-isocenter(2)])
    zlim([0*resolution(3)-isocenter(3),zNumPlanes*resolution(3)-isocenter(3)])
    hold on;
    plot3(x_min,y_min,z_min,'or');
    plot3(x_max,y_max,z_max,'og');
    plot3([x_min x_max],[y_min y_max],[z_min z_max],'w');
    hold off
    grid off
    view(-200,27);
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    title('Y lateral distances')
    colorbar
    box
    
    %%%%%%%%%%%%%%
    subplot(2,3,6)
    z_latDists_vis = zeros(size(ct));
	z_latDists_vis(V) = z_latDists;
    slice(x,y,z,z_latDists_vis,1,1,1) % Draw some volume boundaries
    borders = slice(x,y,z,z_latDists_vis,xslice,yslice,zslice);
    set(borders,'FaceColor','interp',...
        'EdgeColor','none',...
        'DiffuseStrength',.8)
    xlim([0*resolution(1)-isocenter(1),xNumPlanes*resolution(1)-isocenter(1)])
    ylim([0*resolution(2)-isocenter(2),yNumPlanes*resolution(2)-isocenter(2)])
    zlim([0*resolution(3)-isocenter(3),zNumPlanes*resolution(3)-isocenter(3)])
    hold on;
    plot3(x_min,y_min,z_min,'or');
    plot3(x_max,y_max,z_max,'og');
    plot3([x_min x_max],[y_min y_max],[z_min z_max],'w')
    hold off
    grid off
    view(-200,27);
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    title('Z lateral distances')
    colorbar
    box 
    
    figure;
    rad_distances_vis = zeros(size(ct));
	rad_distances_vis(V) = sqrt(rad_distancesSq(ix));
    slice(x,y,z,rad_distances_vis,1,1,1) % Draw some volume boundaries
    borders = slice(x,y,z,rad_distances_vis,xslice,yslice,zslice);
    set(borders,'FaceColor','interp',...
        'EdgeColor','none',...
        'DiffuseStrength',.8)
    xlim([0*resolution(1)-isocenter(1),xNumPlanes*resolution(1)-isocenter(1)])
    ylim([0*resolution(2)-isocenter(2),yNumPlanes*resolution(2)-isocenter(2)])
    zlim([0*resolution(3)-isocenter(3),zNumPlanes*resolution(3)-isocenter(3)])
    hold on;
    plot3(x_min,y_min,z_min,'or');
    plot3(x_max,y_max,z_max,'og');
    plot3([x_min x_max],[y_min y_max],[z_min z_max],'g');
    hold off
    grid off
    view(-200,27);
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    title('Z lateral distances')
    colorbar
    box 
  
    clear coor_x;
    clear coor_y;
    clear coor_z;
    pause;
end