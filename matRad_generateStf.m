function stf = matRad_generateStf(ct,cst,pln,visMode)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad steering information generation
% 
% call
%   stf = matRad_generateStf(ct,cst,pln,visMode)
%
% input
%   ct:         ct cube
%   cst:        matRad cst struct
%   pln:        matRad plan meta information struct
%   visMode:    toggle on/off different visualizations by setting this value to 1,2,3 (optional)
%
% output
%   stf:        matRad steering information struct
%
% References
%   -
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

fprintf('matRad: Generating stf struct... ');

if nargin < 4
    visMode = 0;
end

%
if ~strcmp(pln.radiationMode,'carbon') && sum(strcmp(pln.bioOptimization,{'effect','RBExD'}))>0
    fprintf('\n ********************************************************************************************************* \n');
    fprintf('\n effect based optimization not possible with current setup - physical optimization is carried out instead \n');
    fprintf('\n ********************************************************************************************************* \n');
end

% find all target voxels from cst cell array
V = [];
for i=1:size(cst,1)
    if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6})
        V = [V;cst{i,4}];
    end
end

% Remove double voxels
V = unique(V);

% add margin
addmarginBool = 1;
if addmarginBool
    voi    = zeros(size(ct.cube));
    voi(V) = 1;
    voi    = matRad_addMargin(voi,ct.resolution,ct.resolution,true);
    V      = find(voi>0);
end

% prepare structures necessary for particles
if strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    
    % load base data    
    if strcmp(pln.radiationMode,'protons')
        load protonBaseData;
    elseif  strcmp(pln.radiationMode,'carbon')
        load carbonBaseData;
    end
    
    availableEnergies = [baseData.energy];
    availablePeakPos  = [baseData.peakPos];
    
    clear baseData;
    
end

% Convert linear indices to 3D voxel coordinates
[coordsY_vox, coordsX_vox, coordsZ_vox] = ind2sub(size(ct.cube),V);

% Correct for iso center position. Whit this correction Isocenter is
% (0,0,0) [mm]
coordsX = coordsX_vox*ct.resolution(1) - pln.isoCenter(1);
coordsY = coordsY_vox*ct.resolution(2) - pln.isoCenter(2);
coordsZ = coordsZ_vox*ct.resolution(3) - pln.isoCenter(3);

% Define steering file like struct. Prellocating for speed.
stf = struct;

% loop over all angles
for i = 1:length(pln.gantryAngles)
    
    % Save meta information for treatment plan
    stf(i).gantryAngle   = pln.gantryAngles(i);
    stf(i).couchAngle    = pln.couchAngles(i);
    stf(i).bixelWidth    = pln.bixelWidth;
    stf(i).radiationMode = pln.radiationMode;
    
    % gantry and couch roation matrices according to IEC 61217 standard
    % instead of moving the beam around the patient, we perform an inverse
    % rotation of the patient, i.e. we consider a beam's eye view
    % coordinate system
    
    % Rotation around Z axis (gantry)
    rotMx_XY = [cosd(pln.gantryAngles(i)) -sind(pln.gantryAngles(i)) 0;
                sind(pln.gantryAngles(i))  cosd(pln.gantryAngles(i)) 0;
                                        0                          0 1];
    
    % Rotation around Y axis (Couch movement)
    rotMx_XZ = [ cosd(pln.couchAngles(i)) 0 sind(pln.couchAngles(i));
                                        0 1                         0;
                -sind(pln.couchAngles(i)) 0  cosd(pln.couchAngles(i))];
    
    % rotate target coordinates around Y axis and then around Z axis
    % i.e. 1st couch, 2nd gantry; matrix multiplication not cummutative
    rot_coords = [coordsX coordsY coordsZ]*rotMx_XZ*rotMx_XY;
    
    % project x and z coordinates to isocenter
    coordsAtIsoCenterPlane(:,1) = (rot_coords(:,1)*pln.SAD)./(pln.SAD + rot_coords(:,2));
    coordsAtIsoCenterPlane(:,2) = (rot_coords(:,3)*pln.SAD)./(pln.SAD + rot_coords(:,2));
    
    % Take unique rows values for beamlets positions. Calculate position of
    % central ray for every bixel    
    rayPos = unique(pln.bixelWidth*round([            coordsAtIsoCenterPlane(:,1) ... 
                                          zeros(size(coordsAtIsoCenterPlane,1),1) ...
                                                      coordsAtIsoCenterPlane(:,2)]/pln.bixelWidth),'rows');
    
    % Save the number of rays
    stf(i).numOfRays = size(rayPos,1);
    
    % Save ray and target position in beam eye´s view (bev)
    for j = 1:stf(i).numOfRays
        stf(i).ray(j).rayPos_bev = rayPos(j,:);
        stf(i).ray(j).targetPoint_bev = [2*stf(i).ray(j).rayPos_bev(1) ...
                                                               pln.SAD ...
                                         2*stf(i).ray(j).rayPos_bev(3)];
    end
    
    % source position in bev
    sourcePoint_bev = [0 -pln.SAD 0];
    
    % compute coordinates in lps coordinate system, i.e. rotate beam
    % geometry around fixed patient
    
    % Rotation around Z axis (gantry)
    rotMx_XY_rotated = [ cosd(pln.gantryAngles(i)) sind(pln.gantryAngles(i)) 0;
                        -sind(pln.gantryAngles(i)) cosd(pln.gantryAngles(i)) 0;
                                                 0                         0 1];
    
    % Rotation around Y axis (couch)
    rotMx_XZ_rotated = [ cosd(pln.couchAngles(i)) 0 -sind(pln.couchAngles(i));
                                                0 1                        0;
                         sind(pln.couchAngles(i)) 0 cosd(pln.couchAngles(i))];
    
    % Rotated Source point, first needs to be rotated around gantry, and then
    % couch.
    stf(i).sourcePoint = sourcePoint_bev*rotMx_XY_rotated*rotMx_XZ_rotated;
    
    % Save ray and target position in lps system.
    for j = 1:stf(i).numOfRays
        stf(i).ray(j).rayPos      = stf(i).ray(j).rayPos_bev*rotMx_XY_rotated*rotMx_XZ_rotated;
        stf(i).ray(j).targetPoint = stf(i).ray(j).targetPoint_bev*rotMx_XY_rotated*rotMx_XZ_rotated;
    end
    
    % find appropriate energies for particles
    if strcmp(stf(i).radiationMode,'protons') || strcmp(stf(i).radiationMode,'carbon')
        
        for j = 1:stf(i).numOfRays
            
            % ray tracing necessary to determine depth of the target
            [~,l,rho,~] = matRad_siddonRayTracer(pln.isoCenter,ct.resolution,stf(i).sourcePoint,stf(i).ray(j).targetPoint,{ct.cube,voi});
            
            if sum(rho{2}) > 0 % target hit
                
                % compute radiological depths
                % http://www.ncbi.nlm.nih.gov/pubmed/4000088, eq 14
                radDepths = cumsum(l .* rho{1}); 
                
                % find target entry & exit
                diff_voi    = diff([rho{2}]);
                targetEntry = radDepths(diff_voi == 1);
                targetExit  = radDepths(diff_voi == -1);
                
                if numel(targetEntry) ~= numel(targetExit)
                    error('Inconsistency during ray tracing\n\n');
                end
                
                stf(i).ray(j).energy = [];
                
                % Save energies in stf struct
                for k = 1:numel(targetEntry)
                    stf(i).ray(j).energy = [stf(i).ray(j).energy availableEnergies(availablePeakPos>=targetEntry(k)&availablePeakPos<=targetExit(k))];
                end
                
            else % target not hit
                stf(i).ray(j).energy = [];
            end
            
            % count bixels per ray
            stf(i).numOfBixelsPerRay(j) = numel(stf(i).ray(j).energy);
            
        end
        
    elseif strcmp(stf(i).radiationMode,'photons')
        % set dummy values for photons
        for j = 1:stf(i).numOfRays
            stf(i).ray(j).energy = NaN;
            stf(i).numOfBixelsPerRay(j) = 1;
        end
    else
        error('Error generating stf struct: invalid radiation modality\n\n');
    end
    %%
    
    % save total number of bixels
    stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);
    
    %% visualization
    if visMode > 0
        
        clf;
        % first subplot: visualization in bev
        subplot(1,2,1)
        hold on
        
        % plot rotated target coordinates
        plot3(rot_coords(:,1),rot_coords(:,2),rot_coords(:,3),'r.')
        
        % surface rendering
        if visMode == 2
            
            % generate a 3D rectangular grid centered at isocenter in
            % voxel coordinates
            [X,Y,Z] = meshgrid((1:size(ct.cube,2))-pln.isoCenter(1)/ct.resolution(1), ...
                               (1:size(ct.cube,1))-pln.isoCenter(2)/ct.resolution(2), ...
                               (1:size(ct.cube,3))-pln.isoCenter(3)/ct.resolution(3));
            
            % computes surface
            patSurfCube = 0*ct.cube;
            patSurfCube(unique(cell2mat(cst(:,4)))) = 1;
            [f,v] = isosurface(X,Y,Z,patSurfCube,.5);
            
            % convert isosurface from voxel to [mm]
            v(:,1) = v(:,1)*ct.resolution(1);
            v(:,2) = v(:,2)*ct.resolution(2);
            v(:,3) = v(:,3)*ct.resolution(3);
            
            % rotate surface
            rotated_surface = v*rotMx_XZ*rotMx_XY;
            
            % surface rendering
            surface = patch('Faces',f,'Vertices',rotated_surface);
            set(surface,'FaceColor',[0 0 1],'EdgeColor','none','FaceAlpha',.4);
            lighting gouraud;
        
        end
        
        % plot projection matrix: coordinates at isocenter
        plot3(rayPos(:,1),rayPos(:,2),rayPos(:,3),'k.');
        
        % Plot matrix border of matrix at isocenter
        for j = 1:stf(i).numOfRays
            
            % Compute border for every bixels
            targetPoint_vox_X_1 = stf(i).ray(j).targetPoint_bev(:,1) + pln.bixelWidth;
            targetPoint_vox_Y_1 = stf(i).ray(j).targetPoint_bev(:,2);
            targetPoint_vox_Z_1 = stf(i).ray(j).targetPoint_bev(:,3) + pln.bixelWidth;
            
            targetPoint_vox_X_2 = stf(i).ray(j).targetPoint_bev(:,1) + pln.bixelWidth;
            targetPoint_vox_Y_2 = stf(i).ray(j).targetPoint_bev(:,2);
            targetPoint_vox_Z_2 = stf(i).ray(j).targetPoint_bev(:,3) - pln.bixelWidth;
            
            targetPoint_vox_X_3 = stf(i).ray(j).targetPoint_bev(:,1) - pln.bixelWidth;
            targetPoint_vox_Y_3 = stf(i).ray(j).targetPoint_bev(:,2);
            targetPoint_vox_Z_3 = stf(i).ray(j).targetPoint_bev(:,3) - pln.bixelWidth;
            
            targetPoint_vox_X_4 = stf(i).ray(j).targetPoint_bev(:,1) - pln.bixelWidth;
            targetPoint_vox_Y_4 = stf(i).ray(j).targetPoint_bev(:,2);
            targetPoint_vox_Z_4 = stf(i).ray(j).targetPoint_bev(:,3) + pln.bixelWidth;
            
            % plot
            plot3([sourcePoint_bev(1) targetPoint_vox_X_1],[sourcePoint_bev(2) targetPoint_vox_Y_1],[sourcePoint_bev(3) targetPoint_vox_Z_1],'g')
            plot3([sourcePoint_bev(1) targetPoint_vox_X_2],[sourcePoint_bev(2) targetPoint_vox_Y_2],[sourcePoint_bev(3) targetPoint_vox_Z_2],'g')
            plot3([sourcePoint_bev(1) targetPoint_vox_X_3],[sourcePoint_bev(2) targetPoint_vox_Y_3],[sourcePoint_bev(3) targetPoint_vox_Z_3],'g')
            plot3([sourcePoint_bev(1) targetPoint_vox_X_4],[sourcePoint_bev(2) targetPoint_vox_Y_4],[sourcePoint_bev(3) targetPoint_vox_Z_4],'g')
            
        end
        
        % Plot properties
        daspect([1 1 1]);
        view(0,-90);
        xlabel 'X [mm]'
        ylabel 'Y [mm]'
        zlabel 'Z [mm]'
        title ('Beam''s eye view')
        axis([-300 300 -300 300 -300 300]);
        
        % second subplot: visualization in lps coordinate system
        subplot(1,2,2)
        
        % Plot target coordinates whitout any rotation
        plot3(coordsX,coordsY,coordsZ,'r.')
        hold on;
        
        % Rotated projection matrix at isocenter
        isocenter_plane_coor = rayPos*rotMx_XY_rotated*rotMx_XZ_rotated;
        
        % Plot isocenter plane
        plot3(isocenter_plane_coor(:,1),isocenter_plane_coor(:,2),isocenter_plane_coor(:,3),'y.');
        
        % Plot rotated bixels border.
        for j = 1:stf(i).numOfRays
            % Generate rotated projection target points.
            targetPoint_vox_1_rotated = [stf(i).ray(j).targetPoint_bev(:,1) + pln.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) + pln.bixelWidth]*rotMx_XY_rotated*rotMx_XZ_rotated;
            targetPoint_vox_2_rotated = [stf(i).ray(j).targetPoint_bev(:,1) + pln.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) - pln.bixelWidth]*rotMx_XY_rotated*rotMx_XZ_rotated;
            targetPoint_vox_3_rotated = [stf(i).ray(j).targetPoint_bev(:,1) - pln.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) - pln.bixelWidth]*rotMx_XY_rotated*rotMx_XZ_rotated;
            targetPoint_vox_4_rotated = [stf(i).ray(j).targetPoint_bev(:,1) - pln.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) + pln.bixelWidth]*rotMx_XY_rotated*rotMx_XZ_rotated;
            
            % Plot rotated target points.
            plot3([stf(i).sourcePoint(1) targetPoint_vox_1_rotated(:,1)],[stf(i).sourcePoint(2) targetPoint_vox_1_rotated(:,2)],[stf(i).sourcePoint(3) targetPoint_vox_1_rotated(:,3)],'g')
            plot3([stf(i).sourcePoint(1) targetPoint_vox_2_rotated(:,1)],[stf(i).sourcePoint(2) targetPoint_vox_2_rotated(:,2)],[stf(i).sourcePoint(3) targetPoint_vox_2_rotated(:,3)],'g')
            plot3([stf(i).sourcePoint(1) targetPoint_vox_3_rotated(:,1)],[stf(i).sourcePoint(2) targetPoint_vox_3_rotated(:,2)],[stf(i).sourcePoint(3) targetPoint_vox_3_rotated(:,3)],'g')
            plot3([stf(i).sourcePoint(1) targetPoint_vox_4_rotated(:,1)],[stf(i).sourcePoint(2) targetPoint_vox_4_rotated(:,2)],[stf(i).sourcePoint(3) targetPoint_vox_4_rotated(:,3)],'g')
        end
        
        % surface rendering
        if visMode == 2
            surface = patch('Faces',f,'Vertices',v);
            set(surface,'FaceColor',[0 0 1],'EdgeColor','none','FaceAlpha',.4);
            lighting gouraud;
        end
        
        % labels etc.
        daspect([1 1 1]);
        view(0,-90);
        xlabel 'X [mm]'
        ylabel 'Y [mm]'
        zlabel 'Z [mm]'
        title 'lps coordinate system'
        axis([-300 300 -300 300 -300 300]);
        %pause(1);
    end
    
    % Show progress
    matRad_progress(i,length(pln.gantryAngles));
    
end