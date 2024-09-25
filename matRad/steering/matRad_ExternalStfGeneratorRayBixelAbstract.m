classdef matRad_ExternalStfGeneratorRayBixelAbstract < matRad_StfGeneratorBase
% matRad_PhotonStfGeneratorRayBixelAbstract: Abstract Superclass for 
%   external beam stf generators using the ray-bixel mechanism 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        gantryAngles = 0;
        couchAngles  = 0;
        bixelWidth   = 0;
        isoCenter
        fillEmptyBixels
        centered
    end

    properties (Dependent)
        numOfBeams;
    end

    properties (Access = protected, Hidden)  
        lockAngleUpdate = false;
    end


    methods 
         function this = matRad_ExternalStfGeneratorRayBixelAbstract(pln)
            % Constructs ExternalStfGenerator with or without pln
            if nargin < 1
                pln = [];
            end
            this@matRad_StfGeneratorBase(pln);
         end

         function setDefaults(this)
            % Set default values for ExternalStfGenerator

            this.setDefaults@matRad_StfGeneratorBase();
         end

         function nBeams = get.numOfBeams(this)
            % Return number of beams obtained from angles

            nBeams = numel(this.gantryAngles);
            if nBeams ~= numel(this.couchAngles)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('For some reason, we have a different number of beam and couch angles!');
            end
         end

         function set.gantryAngles(this,angles)
            % Set gantry angles and update couch angles if necessary

            validateattributes(angles,{'numeric'},{'vector','nonempty','nonnan'});            
            oldAngles = this.gantryAngles;
            this.gantryAngles = angles;
            if ~this.lockAngleUpdate
                this.lockAngleUpdate = true;
                if numel(this.gantryAngles) > numel(this.couchAngles)
                    %Append Couch angles with zeros
                    this.couchAngles = [this.couchAngles zeros(1,numel(this.gantryAngles)-numel(this.couchAngles))];                    
                elseif numel(this.couchAngles) > numel(this.gantryAngles)
                    %Try to identify the removed beam angles
                    [removedAngles,ix] = setdiff(oldAngles,this.gantryAngles);
                    
                    nRemovedAngles = numel(this.couchAngles) - numel(this.gantryAngles);
                    
                    if ~isempty(ix) && numel(ix) == nRemovedAngles
                        %Remove corresponding couch angles
                        this.couchAngles(ix) = [];                    
                    else                        
                        this.couchAngles(end-nRemovedAngles+1:end) = [];
                    end
                end
                this.lockAngleUpdate = false;
            end            
         end

         function set.couchAngles(this,angles)
            % Set couch angles and update gantry angles if necessary

            validateattributes(angles,{'numeric'},{'vector','nonempty','nonnan'});            
            oldAngles = this.couchAngles;
            this.couchAngles = angles;
            if ~this.lockAngleUpdate
                this.lockAngleUpdate = true;
                if numel(this.couchAngles) > numel(this.gantryAngles)
                    %Append Gantry angles with zeros
                    this.gantryAngles = [this.gantryAngles zeros(1,numel(this.couchAngles)-numel(this.gantryAngles))];                    
                elseif numel(this.gantryAngles) > numel(this.couchAngles)
                    %Try to identify the removed couch angles
                    [removedAngles,ix] = setdiff(oldAngles,this.couchAngles);

                    nRemovedAngles = numel(this.gantryAngles) - numel(this.couchAngles);

                    if ~isempty(ix) && numel(ix) == nRemovedAngles
                        %Remove corresponding gantry angles
                        this.gantryAngles(ix) = [];                    
                    else                        
                        this.gantryAngles(end-nRemovedAngles+1:end) = [];
                    end                    
                end
                this.lockAngleUpdate = false;
            end            
         end
    end

    methods (Access = protected)
        function initialize(this)
            
            this.initialize@matRad_StfGeneratorBase();

            if this.visMode > 1
                visBool = true;
            else 
                visBool = false;
            end

            if isempty(this.isoCenter)
                this.isoCenter = matRad_getIsoCenter(this.cst,this.ct,visBool);
            end

            if ~isequal(size(this.isoCenter),[this.numOfBeams,3]) && ~size(this.isoCenter,1) ~= 1
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('IsoCenter invalid, creating new one automatically!');
                this.isoCenter = matRad_getIsoCenter(this.cst,this.ct,visBool);
            end
            
            if size(this.isoCenter,1) == 1          
                this.isoCenter = repmat(this.isoCenter,this.numOfBeams,1);
            end
        end

        function rayPos = getRayPositionMatrix(this, beam)
            % Correct for iso center position. Whit this correction Isocenter is
            % (0,0,0) [mm]
            isoCoords = this.voxTargetWorldCoords - beam.isoCenter;
            
            % Get the (active) rotation matrix. We perform a passive/system
            % rotation with row vector coordinates, which would introduce two
            % inversions / transpositions of the matrix, thus no changes to the
            % rotation matrix are necessary           
            rotMat_system_T = matRad_getRotationMatrix(beam.gantryAngle,beam.couchAngle);

            coordsAtIsoCenterPlane =  isoCoords*rotMat_system_T;          
            coordsAtIsoCenterPlane = (coordsAtIsoCenterPlane*beam.SAD)./(beam.SAD+coordsAtIsoCenterPlane(:,2));  
            coordsAtIsoCenterPlane(:,2) = 0;

            rayPos = unique(beam.bixelWidth*round(coordsAtIsoCenterPlane./beam.bixelWidth),'rows');

            % pad ray position array if resolution of target voxel grid not sufficient
            maxCtResolution = max([this.ct.resolution.x this.ct.resolution.y this.ct.resolution.z]);
            if beam.bixelWidth < maxCtResolution
                origRayPos = rayPos;
                for j = -floor(maxCtResolution/beam.bixelWidth):floor(maxCtResolution/beam.bixelWidth)
                    for k = -floor(maxCtResolution/beam.bixelWidth):floor(maxCtResolution/beam.bixelWidth)
                        if abs(j)+abs(k)==0
                            continue;
                        end
                        rayPos = [rayPos; origRayPos(:,1)+j*beam.bixelWidth origRayPos(:,2) origRayPos(:,3)+k*beam.bixelWidth];
                    end
                end
            end

            % remove spaces within rows of bixels for DAO
            %TODO Only Photons
            if this.fillEmptyBixels
                % create single x,y,z vectors
                x = rayPos(:,1);
                y = rayPos(:,2);
                z = rayPos(:,3);
                uniZ = unique(z);
                for j = 1:numel(uniZ)
                    x_loc = x(z == uniZ(j));
                    x_min = min(x_loc);
                    x_max = max(x_loc);
                    x = [x; (x_min:beam.bixelWidth:x_max)'];
                    y = [y; zeros((x_max-x_min)/beam.bixelWidth+1,1)];
                    z = [z; uniZ(j)*ones((x_max-x_min)/beam.bixelWidth+1,1)];
                end

                rayPos = [x,y,z];
            end

             % remove double rays
             rayPos = unique(rayPos,'rows');
        end

        function beam = initBeamData(this,beamIndex)
            beam.gantryAngle     = this.gantryAngles(beamIndex);
            beam.couchAngle      = this.couchAngles(beamIndex);
            beam.isoCenter       = this.isoCenter(beamIndex,:);
            beam.bixelWidth      = this.bixelWidth;
            beam.radiationMode   = this.radiationMode;
            beam.machine         = this.machine.meta.machine;
            beam.SAD             = this.machine.meta.SAD;
        end
            
        
        function stf = generateSourceGeometry(this)
            % Generate basic source geometry (beams & rays) for external beam therapy

            matRad_cfg = MatRad_Config.instance;
   
            % loop over all angles
            for i = 1:length(this.gantryAngles)
                % Save meta information for treatment plan
                beam = this.initBeamData(i);

                beam = this.initRays(beam);                               

                beam = this.setBeamletEnergies(beam);

                if ~isfield(beam.ray,'energy')
                    matRad_cfg.dispError('Error generating stf struct: no suitable energies found. Check if bixelwidth is too large.');
                end                

                beam = this.finalizeBeam(beam); % Post Processing for Ions 

                stf(i) = beam;

                % Show progress
                if matRad_cfg.logLevel > 2
                    matRad_progress(i,length(this.gantryAngles));
                end

                %% visualization
                if this.visMode > 0
                    
                    %center coordinates
                    x = this.ct.x - stf(i).isoCenter(1);
                    y = this.ct.y - stf(i).isoCenter(2);
                    z = this.ct.z - stf(i).isoCenter(3);
                    d = sqrt(sum(([x(1) y(2) z(3)] - [x(end) y(end) z(end)]).^2));
                    limits = [-d/2 d/2 -d/2 d/2 -d/2 d/2];

                    clf;
                    % first subplot: visualization in bev
                    hAxBEV = subplot(1,2,1);
                    set(hAxBEV,'DataAspectRatioMode','manual');
                    hold on;

                    isoCoords = this.voxTargetWorldCoords - stf(i).isoCenter;
            
                    % Get the (active) rotation matrix. We perform a passive/system
                    % rotation with row vector coordinates, which would introduce two
                    % inversions / transpositions of the matrix, thus no changes to the
                    % rotation matrix are necessary           
                    rotMat_system_T = matRad_getRotationMatrix(stf(i).gantryAngle,stf(i).couchAngle);

                    rot_coords = isoCoords*rotMat_system_T;

                    % plot rotated target coordinates
                    plot3(rot_coords(:,1),rot_coords(:,2),rot_coords(:,3),'r.')

                    % surface rendering
                    if this.visMode == 2

                        % computes surface
                        patSurfCube      = 0*this.ct.cube{1};
                        idx              = [this.cst{:,4}];
                        idx              = unique(vertcat(idx{:}));
                        patSurfCube(idx) = 1;

                        [f,v] = isosurface(x,y,z,patSurfCube,.5);

                        vRot = v*rotMat_system_T;
    
                        % rotate surface
                        rotated_surface = v*rotMat_system_T;

                        % surface rendering
                        surface = patch('Faces',f,'Vertices',rotated_surface);
                        set(surface,'FaceColor',[0 0 1],'EdgeColor','none','FaceAlpha',.4);
                        lighting gouraud;

                    end

                    % plot projection matrix: coordinates at isocenter
                    plot3(rayPos_bev(:,1),rayPos_bev(:,2),rayPos_bev(:,3),'k.');

                    % Plot matrix border of matrix at isocenter
                    for j = 1:stf(i).numOfRays

                        % Compute border for every bixels
                        targetPoint_vox_X_1 = stf(i).ray(j).targetPoint_bev(:,1) + stf(i).bixelWidth;
                        targetPoint_vox_Y_1 = stf(i).ray(j).targetPoint_bev(:,2);
                        targetPoint_vox_Z_1 = stf(i).ray(j).targetPoint_bev(:,3) + stf(i).bixelWidth;

                        targetPoint_vox_X_2 = stf(i).ray(j).targetPoint_bev(:,1) + stf(i).bixelWidth;
                        targetPoint_vox_Y_2 = stf(i).ray(j).targetPoint_bev(:,2);
                        targetPoint_vox_Z_2 = stf(i).ray(j).targetPoint_bev(:,3) - stf(i).bixelWidth;

                        targetPoint_vox_X_3 = stf(i).ray(j).targetPoint_bev(:,1) - stf(i).bixelWidth;
                        targetPoint_vox_Y_3 = stf(i).ray(j).targetPoint_bev(:,2);
                        targetPoint_vox_Z_3 = stf(i).ray(j).targetPoint_bev(:,3) - stf(i).bixelWidth;

                        targetPoint_vox_X_4 = stf(i).ray(j).targetPoint_bev(:,1) - stf(i).bixelWidth;
                        targetPoint_vox_Y_4 = stf(i).ray(j).targetPoint_bev(:,2);
                        targetPoint_vox_Z_4 = stf(i).ray(j).targetPoint_bev(:,3) + stf(i).bixelWidth;

                        % plot
                        plot3([stf(i).sourcePoint_bev(1) targetPoint_vox_X_1],[stf(i).sourcePoint_bev(2) targetPoint_vox_Y_1],[stf(i).sourcePoint_bev(3) targetPoint_vox_Z_1],'g')
                        plot3([stf(i).sourcePoint_bev(1) targetPoint_vox_X_2],[stf(i).sourcePoint_bev(2) targetPoint_vox_Y_2],[stf(i).sourcePoint_bev(3) targetPoint_vox_Z_2],'g')
                        plot3([stf(i).sourcePoint_bev(1) targetPoint_vox_X_3],[stf(i).sourcePoint_bev(2) targetPoint_vox_Y_3],[stf(i).sourcePoint_bev(3) targetPoint_vox_Z_3],'g')
                        plot3([stf(i).sourcePoint_bev(1) targetPoint_vox_X_4],[stf(i).sourcePoint_bev(2) targetPoint_vox_Y_4],[stf(i).sourcePoint_bev(3) targetPoint_vox_Z_4],'g')

                    end

                    % Plot properties
                    view(hAxBEV,0,-90);
                    xlabel(hAxBEV,'X [mm]');
                    ylabel(hAxBEV,'Y [mm]');
                    zlabel(hAxBEV,'Z [mm]');
                    title(hAxBEV,'Beam''s eye view');

                    axis(hAxBEV,limits);

                    % second subplot: visualization in lps coordinate system
                    hAxLPS = subplot(1,2,2);

                    % Plot target coordinates whitout any rotation
                    plot3(isoCoords(:,1),isoCoords(:,2),isoCoords(:,3),'r.')
                    hold on;

                    % Rotated projection matrix at isocenter
                    isocenter_plane_coor = rayPos_bev*rotMat_vectors_T;

                    % Plot isocenter plane
                    plot3(isocenter_plane_coor(:,1),isocenter_plane_coor(:,2),isocenter_plane_coor(:,3),'y.');

                    % Plot rotated bixels border.
                    for j = 1:stf(i).numOfRays
                        % Generate rotated projection target points.
                        targetPoint_vox_1_rotated = [stf(i).ray(j).targetPoint_bev(:,1) + stf(i).bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) + stf(i).bixelWidth]*rotMat_vectors_T;
                        targetPoint_vox_2_rotated = [stf(i).ray(j).targetPoint_bev(:,1) + stf(i).bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) - stf(i).bixelWidth]*rotMat_vectors_T;
                        targetPoint_vox_3_rotated = [stf(i).ray(j).targetPoint_bev(:,1) - stf(i).bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) - stf(i).bixelWidth]*rotMat_vectors_T;
                        targetPoint_vox_4_rotated = [stf(i).ray(j).targetPoint_bev(:,1) - stf(i).bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) + stf(i).bixelWidth]*rotMat_vectors_T;

                        % Plot rotated target points.
                        plot3([stf(i).sourcePoint(1) targetPoint_vox_1_rotated(:,1)],[stf(i).sourcePoint(2) targetPoint_vox_1_rotated(:,2)],[stf(i).sourcePoint(3) targetPoint_vox_1_rotated(:,3)],'g')
                        plot3([stf(i).sourcePoint(1) targetPoint_vox_2_rotated(:,1)],[stf(i).sourcePoint(2) targetPoint_vox_2_rotated(:,2)],[stf(i).sourcePoint(3) targetPoint_vox_2_rotated(:,3)],'g')
                        plot3([stf(i).sourcePoint(1) targetPoint_vox_3_rotated(:,1)],[stf(i).sourcePoint(2) targetPoint_vox_3_rotated(:,2)],[stf(i).sourcePoint(3) targetPoint_vox_3_rotated(:,3)],'g')
                        plot3([stf(i).sourcePoint(1) targetPoint_vox_4_rotated(:,1)],[stf(i).sourcePoint(2) targetPoint_vox_4_rotated(:,2)],[stf(i).sourcePoint(3) targetPoint_vox_4_rotated(:,3)],'g')
                    end

                    % surface rendering
                    if this.visMode == 2
                        surface = patch('Faces',f,'Vertices',v);
                        set(surface,'FaceColor',[0 0 1],'EdgeColor','none','FaceAlpha',.4);
                        lighting gouraud;
                    end

                    % labels etc.
                    daspect([1 1 1]);
                    view(0,-90);
                    xlabel(hAxBEV,'X [mm]');
                    ylabel(hAxBEV,'Y [mm]');
                    zlabel(hAxBEV,'Z [mm]');
                    title('lps coordinate system');
                    axis(limits);
                    drawnow();
                    pause(1);
                end
                

                % Show progress
                if matRad_cfg.logLevel > 2
                    matRad_progress(i,length(this.gantryAngles));
                end
            end
        end

        function beam = initRays(this,beam)
            % source position in bev
            beam.sourcePoint_bev = [0 -beam.SAD 0];

            rayPos_bev = this.getRayPositionMatrix(beam);

            % Save the number of rays
            beam.numOfRays = size(rayPos_bev,1);

            % Save ray and target position in beam eye's view (bev)
            for j = 1:beam.numOfRays
                beam.ray(j).rayPos_bev = rayPos_bev(j,:);
                beam.ray(j).targetPoint_bev = [2*beam.ray(j).rayPos_bev(1) ...
                    beam.SAD ...
                    2*beam.ray(j).rayPos_bev(3)];
            end

            % get (active) rotation matrix
            % transpose matrix because we are working with row vectors
            rotMat_vectors_T = transpose(matRad_getRotationMatrix(beam.gantryAngle,beam.couchAngle));

            %Source point in LPS system
            beam.sourcePoint = beam.sourcePoint_bev*rotMat_vectors_T;
               
            % Save ray and target position in lps system.
            for j = 1:beam.numOfRays
                beam.ray(j).rayPos      = beam.ray(j).rayPos_bev*rotMat_vectors_T;
                beam.ray(j).targetPoint = beam.ray(j).targetPoint_bev*rotMat_vectors_T;
            end                        
        end

        function beam = setBeamletEnergies(this,beam)
            % Abstract method to set source energy on beam. Implemented in Subclasses

            throw(MException('MATLAB:class:AbstractMember','This method is not implemented in the base class'));
        end


        function beam = finalizeBeam(this,beam)
            %Finalize meta data of beam            

            matRad_cfg = MatRad_Config.instance();

            numOfRays = numel(beam.ray);
            if ~isfield(beam,'numOfRays')
                beam.numOfRays = numOfRays;
            end
            if ~isequal(numOfRays,beam.numOfRays)
                matRad_cfg.dispError('Validation of number of rays failed!');
            end

            % Calculate the number of bixels per ray
            numOfBixelsPerRay = arrayfun(@(r) numel(r.energy),beam.ray);     
            if ~isfield(beam,'numOfBixelsPerRay')
                beam.numOfBixelsPerRay = numOfBixelsPerRay;
            end
            if ~isequal(numOfBixelsPerRay,beam.numOfBixelsPerRay)
                matRad_cfg.dispError('Validation of number of bixels per ray failed!');
            end
            
            % Calculate the total number of bixels            
            totalNumOfBixels = sum(beam.numOfBixelsPerRay);
            if ~isfield(beam,'totalNumOfBixels')
                beam.totalNumOfBixels = totalNumOfBixels;
            end
            if ~isequal(totalNumOfBixels,beam.totalNumOfBixels)
                matRad_cfg.dispError('Validation of total number of bixels failed!');
            end
            
        end
    end

    methods (Static)
        function [available,msg] = isAvailable(pln,machine)
            % see superclass for information            
                   
            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            %checkBasic
            available = isfield(machine,'meta') && isfield(machine,'data');

            available = available && isfield(machine.meta,'machine');

            available = available && isfield(machine.meta,'SAD') && isscalar(machine.meta.SAD);
    
            if ~available
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';                
            else
                msg = [];
            end
        end
    end
end

