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


fprintf('matRad: Generating stf struct... ');

if nargin < 4
    visMode = 0;
end

if numel(pln.propStf.gantryAngles) ~= numel(pln.propStf.couchAngles)
    error('Inconsistent number of gantry and couch angles.');
end

if pln.propStf.bixelWidth < 0 || ~isfinite(pln.propStf.bixelWidth)
   error('bixel width (spot distance) needs to be a real number [mm] larger than zero.');
end

% find all target voxels from cst cell array
V = [];
for i=1:size(cst,1)
    if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6})
        V = [V;vertcat(cst{i,4}{:})];
    end
end

% Remove double voxels
V = unique(V);

% generate voi cube for targets
voiTarget    = zeros(ct.cubeDim);
voiTarget(V) = 1;
    
% add margin
addmarginBool = 1;
if addmarginBool
    voiTarget = matRad_addMargin(voiTarget,cst,ct.resolution,ct.resolution,true);
    V   = find(voiTarget>0);
end

% throw error message if no target is found
if isempty(V)
    error('Could not find target.');
end

% prepare structures necessary for particles
fileName = [pln.radiationMode '_' pln.machine];
try
   load([fileparts(mfilename('fullpath')) filesep fileName]);
   SAD = machine.meta.SAD;
catch
   error(['Could not find the following machine file: ' fileName ]); 
end

if strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
      
    availableEnergies = [machine.data.energy];
    availablePeakPos  = [machine.data.peakPos] + [machine.data.offset];
    
    if sum(availablePeakPos<0)>0
       error('at least one available peak position is negative - inconsistent machine file') 
    end
    %clear machine;
end

% calculate rED or rSP from HU
if ~isdeployed
   addpath(['dicomImport'])
   addpath(['dicomImport' filesep 'hlutLibrary'])
end
ct = matRad_calcWaterEqD(ct, pln);

% Convert linear indices to 3D voxel coordinates
[coordsY_vox, coordsX_vox, coordsZ_vox] = ind2sub(ct.cubeDim,V);

% Define steering file like struct. Prellocating for speed.
stf = struct;

% loop over all angles
for i = 1:length(pln.propStf.gantryAngles)
    
    % Correct for iso center position. Whit this correction Isocenter is
    % (0,0,0) [mm]
    coordsX = coordsX_vox*ct.resolution.x - pln.propStf.isoCenter(i,1);
    coordsY = coordsY_vox*ct.resolution.y - pln.propStf.isoCenter(i,2);
    coordsZ = coordsZ_vox*ct.resolution.z - pln.propStf.isoCenter(i,3);

    % Save meta information for treatment plan
    stf(i).gantryAngle   = pln.propStf.gantryAngles(i);
    stf(i).couchAngle    = pln.propStf.couchAngles(i);
    stf(i).bixelWidth    = pln.propStf.bixelWidth;
    stf(i).radiationMode = pln.radiationMode;
    stf(i).SAD           = SAD;
    stf(i).isoCenter     = pln.propStf.isoCenter(i,:);
        
    % Get the (active) rotation matrix. We perform a passive/system 
    % rotation with row vector coordinates, which would introduce two 
    % inversions / transpositions of the matrix, thus no changes to the
    % rotation matrix are necessary
    rotMat_system_T = matRad_getRotationMatrix(pln.propStf.gantryAngles(i),pln.propStf.couchAngles(i));
    
    rot_coords = [coordsX coordsY coordsZ]*rotMat_system_T;
    
    % project x and z coordinates to isocenter
    coordsAtIsoCenterPlane(:,1) = (rot_coords(:,1)*SAD)./(SAD + rot_coords(:,2));
    coordsAtIsoCenterPlane(:,2) = (rot_coords(:,3)*SAD)./(SAD + rot_coords(:,2));
    
    % Take unique rows values for beamlets positions. Calculate position of
    % central ray for every bixel    
    rayPos = unique(pln.propStf.bixelWidth*round([           coordsAtIsoCenterPlane(:,1) ... 
                                                  zeros(size(coordsAtIsoCenterPlane,1),1) ...
                                                             coordsAtIsoCenterPlane(:,2)]/pln.propStf.bixelWidth),'rows');
                                                  
    % pad ray position array if resolution of target voxel grid not sufficient
    maxCtResolution = max([ct.resolution.x ct.resolution.y ct.resolution.z]);
    if pln.propStf.bixelWidth < maxCtResolution
        origRayPos = rayPos;
        for j = -floor(maxCtResolution/pln.propStf.bixelWidth):floor(maxCtResolution/pln.propStf.bixelWidth)
            for k = -floor(maxCtResolution/pln.propStf.bixelWidth):floor(maxCtResolution/pln.propStf.bixelWidth)
                if abs(j)+abs(k)==0
                    continue;
                end                
                rayPos = [rayPos; origRayPos(:,1)+j*pln.propStf.bixelWidth origRayPos(:,2) origRayPos(:,3)+k*pln.propStf.bixelWidth];
            end
        end
     end

     % remove spaces within rows of bixels for DAO
     if pln.propOpt.runDAO
         % create single x,y,z vectors
         x = rayPos(:,1);
         y = rayPos(:,2);
         z = rayPos(:,3);
         uniZ = unique(z);
         for j = 1:numel(uniZ)
             x_loc = x(z == uniZ(j));
             x_min = min(x_loc);
             x_max = max(x_loc);
             x = [x; [x_min:pln.propStf.bixelWidth:x_max]'];
             y = [y; zeros((x_max-x_min)/pln.propStf.bixelWidth+1,1)];
             z = [z; uniZ(j)*ones((x_max-x_min)/pln.propStf.bixelWidth+1,1)];             
         end
         
         rayPos = [x,y,z];
     end
    
    % remove double rays
    rayPos = unique(rayPos,'rows');
    
    % Save the number of rays
    stf(i).numOfRays = size(rayPos,1);
    
    % Save ray and target position in beam eye's view (bev)
    for j = 1:stf(i).numOfRays
        stf(i).ray(j).rayPos_bev = rayPos(j,:);
        stf(i).ray(j).targetPoint_bev = [2*stf(i).ray(j).rayPos_bev(1) ...
                                                               SAD ...
                                         2*stf(i).ray(j).rayPos_bev(3)];
    end
    
    % source position in bev
    stf(i).sourcePoint_bev = [0 -SAD 0];
    
    % get (active) rotation matrix 
    % transpose matrix because we are working with row vectors
    rotMat_vectors_T = transpose(matRad_getRotationMatrix(pln.propStf.gantryAngles(i),pln.propStf.couchAngles(i)));
    
    
    stf(i).sourcePoint = stf(i).sourcePoint_bev*rotMat_vectors_T;
    
    % Save ray and target position in lps system.
    for j = 1:stf(i).numOfRays
        stf(i).ray(j).rayPos      = stf(i).ray(j).rayPos_bev*rotMat_vectors_T;
        stf(i).ray(j).targetPoint = stf(i).ray(j).targetPoint_bev*rotMat_vectors_T;
        if strcmp(pln.radiationMode,'photons') 
            stf(i).ray(j).rayCorners_SCD = (repmat([0, machine.meta.SCD - SAD, 0],4,1)+ (machine.meta.SCD/SAD) * ...
                                                             [rayPos(j,:) + [+stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                                                              rayPos(j,:) + [-stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                                                              rayPos(j,:) + [-stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2];...
                                                              rayPos(j,:) + [+stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2]])*rotMat_vectors_T;
        end
    end
    
    % loop over all rays to determine meta information for each ray    
    stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);
    
    for j = stf(i).numOfRays:-1:1

        % ray tracing necessary to determine depth of the target
        [~,l,rho,~,~] = matRad_siddonRayTracer(stf(i).isoCenter, ...
                             ct.resolution, ...
                             stf(i).sourcePoint, ...
                             stf(i).ray(j).targetPoint, ...
                             [{ct.cube{1}} {voiTarget}]);

        % find appropriate energies for particles
       if strcmp(stf(i).radiationMode,'protons') || strcmp(stf(i).radiationMode,'carbon')

           % target hit
           if sum(rho{2}) > 0 

                % compute radiological depths
                % http://www.ncbi.nlm.nih.gov/pubmed/4000088, eq 14
                radDepths = cumsum(l .* rho{1}); 

                % find target entry & exit
                diff_voi    = diff([rho{2}]);
                targetEntry = radDepths(diff_voi == 1);
                targetExit  = radDepths(diff_voi == -1);

                if numel(targetEntry) ~= numel(targetExit)
                    error('Inconsistency during ray tracing.');
                end

                stf(i).ray(j).energy = [];

                % Save energies in stf struct
                for k = 1:numel(targetEntry)
                    stf(i).ray(j).energy = [stf(i).ray(j).energy availableEnergies(availablePeakPos>=targetEntry(k)&availablePeakPos<=targetExit(k))];
                end
  
                % book keeping & calculate focus index
                stf(i).numOfBixelsPerRay(j) = numel([stf(i).ray(j).energy]);
                currentMinimumFWHM = matRad_interp1(machine.meta.LUT_bxWidthminFWHM(1,:)',...
                                             machine.meta.LUT_bxWidthminFWHM(2,:)',...
                                             pln.propStf.bixelWidth);
                focusIx  =  ones(stf(i).numOfBixelsPerRay(j),1);
                [~, vEnergyIx] = min(abs(bsxfun(@minus,[machine.data.energy]',...
                                repmat(stf(i).ray(j).energy,length([machine.data]),1))));

                % get for each spot the focus index
                for k = 1:stf(i).numOfBixelsPerRay(j)                    
                    focusIx(k) = find(machine.data(vEnergyIx(k)).initFocus.SisFWHMAtIso > currentMinimumFWHM,1,'first');
                end

                stf(i).ray(j).focusIx = focusIx';
                 
            else % target not hit
                stf(i).ray(j)               = [];
                stf(i).numOfBixelsPerRay(j) = [];
           end
           
       elseif strcmp(stf(i).radiationMode,'photons')
           
         % book keeping for photons
         stf(i).ray(j).energy = machine.data.energy;
         
       else
          error('Error generating stf struct: invalid radiation modality.');
       end
       
    end
    
    % store total number of rays for beam-i
    stf(i).numOfRays = size(stf(i).ray,2);
     
    % post processing for particle remove energy slices
    if strcmp(stf(i).radiationMode,'protons') || strcmp(stf(i).radiationMode,'carbon')
        
        % get minimum energy per field
        minEnergy = min([stf(i).ray.energy]);
        maxEnergy = max([stf(i).ray.energy]);
        
        % get corresponding peak position
        availableEnergies = [machine.data.energy];
        minPeakPos  = machine.data(minEnergy == availableEnergies).peakPos;
        maxPeakPos  = machine.data(maxEnergy == availableEnergies).peakPos;
        
        % find set of energyies with adequate spacing
        if ~isfield(pln.propStf, 'longitudinalSpotSpacing')
            if strcmp(machine.meta.machine,'Generic')
                longitudinalSpotSpacing = 1.5; % enforce all entries to be used
            else
                longitudinalSpotSpacing = 3;   % default value for all other treatment machines
            end
        else
            longitudinalSpotSpacing = pln.propStf.longitudinalSpotSpacing;
        end
        
        stf(i).longitudinalSpotSpacing = longitudinalSpotSpacing;
        
        tolerance              = longitudinalSpotSpacing/10;
        availablePeakPos       = [machine.data.peakPos];
        
        useEnergyBool = availablePeakPos >= minPeakPos & availablePeakPos <= maxPeakPos;
        
        ixCurr = find(useEnergyBool,1,'first');
        ixRun  = ixCurr + 1;
        ixEnd  = find(useEnergyBool,1,'last');

        while ixRun <= ixEnd
            if abs(availablePeakPos(ixRun)-availablePeakPos(ixCurr)) < ...
                                    longitudinalSpotSpacing - tolerance
                useEnergyBool(ixRun) = 0;
            else
                ixCurr = ixRun;
            end
            ixRun = ixRun + 1;
        end
        
        for j = stf(i).numOfRays:-1:1
            for k = stf(i).numOfBixelsPerRay(j):-1:1
                maskEnergy = stf(i).ray(j).energy(k) == availableEnergies;
                if ~useEnergyBool(maskEnergy)
                    stf(i).ray(j).energy(k)     = [];
                    stf(i).ray(j).focusIx(k)    = [];
                    stf(i).numOfBixelsPerRay(j) = stf(i).numOfBixelsPerRay(j) - 1;
                end
            end
            if isempty(stf(i).ray(j).energy)
                stf(i).ray(j) = [];
                stf(i).numOfBixelsPerRay(j) = [];
                stf(i).numOfRays = stf(i).numOfRays - 1;
            end
        end
        
    end
    
    % save total number of bixels
    stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);
    
    % Show progress
    matRad_progress(i,length(pln.propStf.gantryAngles));

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
            [X,Y,Z] = meshgrid((1:ct.cubeDim(2))-stf(i).isoCenter(1)/ct.resolution.x, ...
                               (1:ct.cubeDim(1))-stf(i).isoCenter(2)/ct.resolution.y, ...
                               (1:ct.cubeDim(3))-stf(i).isoCenter(3)/ct.resolution.z);
            
            % computes surface
            patSurfCube      = 0*ct.cube{1};
            idx              = [cst{:,4}];
            idx              = unique(vertcat(idx{:}));
            patSurfCube(idx) = 1;
            
            [f,v] = isosurface(X,Y,Z,patSurfCube,.5);
            
            % convert isosurface from voxel to [mm]
            v(:,1) = v(:,1)*ct.resolution.x;
            v(:,2) = v(:,2)*ct.resolution.y;
            v(:,3) = v(:,3)*ct.resolution.z;
            
            % rotate surface
            rotated_surface = v*rotMat_system_T;
            
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
            targetPoint_vox_X_1 = stf(i).ray(j).targetPoint_bev(:,1) + pln.propStf.bixelWidth;
            targetPoint_vox_Y_1 = stf(i).ray(j).targetPoint_bev(:,2);
            targetPoint_vox_Z_1 = stf(i).ray(j).targetPoint_bev(:,3) + pln.propStf.bixelWidth;
            
            targetPoint_vox_X_2 = stf(i).ray(j).targetPoint_bev(:,1) + pln.propStf.bixelWidth;
            targetPoint_vox_Y_2 = stf(i).ray(j).targetPoint_bev(:,2);
            targetPoint_vox_Z_2 = stf(i).ray(j).targetPoint_bev(:,3) - pln.propStf.bixelWidth;
            
            targetPoint_vox_X_3 = stf(i).ray(j).targetPoint_bev(:,1) - pln.propStf.bixelWidth;
            targetPoint_vox_Y_3 = stf(i).ray(j).targetPoint_bev(:,2);
            targetPoint_vox_Z_3 = stf(i).ray(j).targetPoint_bev(:,3) - pln.propStf.bixelWidth;
            
            targetPoint_vox_X_4 = stf(i).ray(j).targetPoint_bev(:,1) - pln.propStf.bixelWidth;
            targetPoint_vox_Y_4 = stf(i).ray(j).targetPoint_bev(:,2);
            targetPoint_vox_Z_4 = stf(i).ray(j).targetPoint_bev(:,3) + pln.propStf.bixelWidth;
            
            % plot
            plot3([stf(i).sourcePoint_bev(1) targetPoint_vox_X_1],[stf(i).sourcePoint_bev(2) targetPoint_vox_Y_1],[stf(i).sourcePoint_bev(3) targetPoint_vox_Z_1],'g')
            plot3([stf(i).sourcePoint_bev(1) targetPoint_vox_X_2],[stf(i).sourcePoint_bev(2) targetPoint_vox_Y_2],[stf(i).sourcePoint_bev(3) targetPoint_vox_Z_2],'g')
            plot3([stf(i).sourcePoint_bev(1) targetPoint_vox_X_3],[stf(i).sourcePoint_bev(2) targetPoint_vox_Y_3],[stf(i).sourcePoint_bev(3) targetPoint_vox_Z_3],'g')
            plot3([stf(i).sourcePoint_bev(1) targetPoint_vox_X_4],[stf(i).sourcePoint_bev(2) targetPoint_vox_Y_4],[stf(i).sourcePoint_bev(3) targetPoint_vox_Z_4],'g')
            
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
        isocenter_plane_coor = rayPos*rotMat_vectors_T;
        
        % Plot isocenter plane
        plot3(isocenter_plane_coor(:,1),isocenter_plane_coor(:,2),isocenter_plane_coor(:,3),'y.');
        
        % Plot rotated bixels border.
        for j = 1:stf(i).numOfRays
            % Generate rotated projection target points.
            targetPoint_vox_1_rotated = [stf(i).ray(j).targetPoint_bev(:,1) + pln.propStf.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) + pln.propStf.bixelWidth]*rotMat_vectors_T;
            targetPoint_vox_2_rotated = [stf(i).ray(j).targetPoint_bev(:,1) + pln.propStf.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) - pln.propStf.bixelWidth]*rotMat_vectors_T;
            targetPoint_vox_3_rotated = [stf(i).ray(j).targetPoint_bev(:,1) - pln.propStf.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) - pln.propStf.bixelWidth]*rotMat_vectors_T;
            targetPoint_vox_4_rotated = [stf(i).ray(j).targetPoint_bev(:,1) - pln.propStf.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) + pln.propStf.bixelWidth]*rotMat_vectors_T;
            
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
    
    % include rangeshifter data if not yet available 
    if strcmp(pln.radiationMode, 'protons') || strcmp(pln.radiationMode, 'carbon')
        for j = 1:stf(i).numOfRays
            for k = 1:numel(stf(i).ray(j).energy)
                stf(i).ray(j).rangeShifter(k).ID = 0;
                stf(i).ray(j).rangeShifter(k).eqThickness = 0;
                stf(i).ray(j).rangeShifter(k).sourceRashiDistance = 0;
            end
        end
    end
        
end    

end
