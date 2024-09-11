classdef matRad_ExternalStfGeneratorIMRT < matRad_StfGeneratorBase

    
    properties (Access = protected)
        rayPos
        ctEntryPoint        
    end

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
         function this = matRad_ExternalStfGeneratorIMRT(pln)
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
        function initializePatientGeometry(this)
            % Initialize patient geometry for external beam therapy geometry

            matRad_cfg = MatRad_Config.instance();

            if this.visMode > 1
                visBool = true;
            else 
                visBool = false;
            end

            if isempty(this.isoCenter)
                this.isoCenter = matRad_getIsoCenter(this.cst,this.ct,visBool);
            end

            if ~isequal(size(this.isoCenter),[this.numOfBeams,3]) && ~size(this.isoCenter,1) ~= 1
                matRad_cfg.dispWarning('IsoCenter invalid, creating new one automatically!');
                this.isoCenter = matRad_getIsoCenter(this.cst,this.ct,visBool);
            end
            
            if size(this.isoCenter,1) ~= 1          
                this.isoCenter = repmat(this.isoCenter,this.numOfBeams,1);
            end

            initializePatientGeometry@matRad_StfGeneratorBase(this);
        end
        
        function stf = generateSourceGeometry(this)
            % Generate basic source geometry (beams & rays) for external beam therapy

            matRad_cfg = MatRad_Config.instance;

            % prepare structures necessary for particles
            SAD = this.machine.meta.SAD;

            % Define steering file like struct. Prellocating for speed.
            stf = struct;
    
            % loop over all angles
            for i = 1:length(this.gantryAngles)
                % Save meta information for treatment plan
                stf(i).gantryAngle   = this.gantryAngles(i);
                stf(i).couchAngle    = this.couchAngles(i);
                stf(i).bixelWidth    = this.bixelWidth;
                stf(i).radiationMode = this.radiationMode;
                stf(i).machine       = this.machine.meta.machine;
                stf(i).SAD           = SAD;
                stf(i).isoCenter     = this.isoCenter(i,:);

                % Correct for iso center position. Whit this correction Isocenter is
                % (0,0,0) [mm]
                isoCoords = this.voxTargetWorldCoords - stf(i).isoCenter;

                % Get the (active) rotation matrix. We perform a passive/system
                % rotation with row vector coordinates, which would introduce two
                % inversions / transpositions of the matrix, thus no changes to the
                % rotation matrix are necessary
                rotMat_system_T = matRad_getRotationMatrix(stf(i).gantryAngle,stf(i).couchAngle);

                rot_coords =  isoCoords*rotMat_system_T;

                % project x and z coordinates to isocenter
                coordsAtIsoCenterPlane(:,1) = (rot_coords(:,1)*SAD)./(SAD + rot_coords(:,2));
                coordsAtIsoCenterPlane(:,2) = (rot_coords(:,3)*SAD)./(SAD + rot_coords(:,2));

                % Take unique rows values for beamlets positions. Calculate position of
                % central ray for every bixel
                this.rayPos = unique(this.bixelWidth*round([coordsAtIsoCenterPlane(:,1) ...
                    zeros(size(coordsAtIsoCenterPlane,1),1) ...
                    coordsAtIsoCenterPlane(:,2)]/this.bixelWidth),'rows');

                % pad ray position array if resolution of target voxel grid not sufficient
                maxCtResolution = max([this.ct.resolution.x this.ct.resolution.y this.ct.resolution.z]);
                if this.bixelWidth < maxCtResolution
                    origRayPos = this.rayPos;
                    for j = -floor(maxCtResolution/this.bixelWidth):floor(maxCtResolution/this.bixelWidth)
                        for k = -floor(maxCtResolution/this.bixelWidth):floor(maxCtResolution/this.bixelWidth)
                            if abs(j)+abs(k)==0
                                continue;
                            end
                            this.rayPos = [this.rayPos; origRayPos(:,1)+j*this.bixelWidth origRayPos(:,2) origRayPos(:,3)+k*this.bixelWidth];
                        end
                    end
                end

                % remove spaces within rows of bixels for DAO
                %TODO Only Photons
                if this.fillEmptyBixels
                    % create single x,y,z vectors
                    x = this.rayPos(:,1);
                    y = this.rayPos(:,2);
                    z = this.rayPos(:,3);
                    uniZ = unique(z);
                    for j = 1:numel(uniZ)
                        x_loc = x(z == uniZ(j));
                        x_min = min(x_loc);
                        x_max = max(x_loc);
                        x = [x; (x_min:this.bixelWidth:x_max)'];
                        y = [y; zeros((x_max-x_min)/this.bixelWidth+1,1)];
                        z = [z; uniZ(j)*ones((x_max-x_min)/this.bixelWidth+1,1)];
                    end

                    this.rayPos = [x,y,z];
                end

                % remove double rays
                this.rayPos = unique(this.rayPos,'rows');

                % Save the number of rays
                stf(i).numOfRays = size(this.rayPos,1);

                % Save ray and target position in beam eye's view (bev)
                for j = 1:stf(i).numOfRays
                    stf(i).ray(j).rayPos_bev = this.rayPos(j,:);
                    stf(i).ray(j).targetPoint_bev = [2*stf(i).ray(j).rayPos_bev(1) ...
                        SAD ...
                        2*stf(i).ray(j).rayPos_bev(3)];
                end

                % source position in bev
                stf(i).sourcePoint_bev = [0 -SAD 0];

                % get (active) rotation matrix
                % transpose matrix because we are working with row vectors
                rotMat_vectors_T = transpose(matRad_getRotationMatrix(this.gantryAngles(i),this.couchAngles(i)));


                stf(i).sourcePoint = stf(i).sourcePoint_bev*rotMat_vectors_T;

                stf(i) = this.initializeRayTargetPosition(stf(i),rotMat_vectors_T,SAD);
                
                % loop over all rays to determine meta information for each ray
                stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);

                stf(i) = this.setSourceEnergyOnBeam(stf(i));

                if ~isfield(stf(i).ray,'energy')
                    matRad_cfg.dispError('Error generating stf struct: no suitable energies found. Check if bixelwidth is too large.');
                end
                % store total number of rays for beam-i
                stf(i).numOfRays = size(stf(i).ray,2);

                stf(i).longitudinalSpotSpacing = []; % ensures the stf(i) structs have the same structures when the next line comes
                stf(i) = this.finalizeBeam(stf(i)); % Post Processing for Ions 

                % save total number of bixels
                stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);

                % Show progress
                if matRad_cfg.logLevel > 2
                    matRad_progress(i,length(this.gantryAngles));
                end

                %% visualization
                if this.visMode > 0

                    clf;
                    % first subplot: visualization in bev
                    subplot(1,2,1)
                    hold on

                    % plot rotated target coordinates
                    plot3(rot_coords(:,1),rot_coords(:,2),rot_coords(:,3),'r.')

                    % surface rendering
                    if this.visMode == 2

                        % generate a 3D rectangular grid centered at isocenter in
                        % voxel coordinates
                        [X,Y,Z] = meshgrid((1:this.ct.cubeDim(2))-stf(i).isoCenter(1)/this.ct.resolution.x, ...
                            (1:this.ct.cubeDim(1))-stf(i).isoCenter(2)/this.ct.resolution.y, ...
                            (1:this.ct.cubeDim(3))-stf(i).isoCenter(3)/this.ct.resolution.z);

                        % computes surface
                        patSurfCube      = 0*this.ct.cube{1};
                        idx              = [this.cst{:,4}];
                        idx              = unique(vertcat(idx{:}));
                        patSurfCube(idx) = 1;

                        [f,v] = isosurface(X,Y,Z,patSurfCube,.5);

                        % convert isosurface from voxel to [mm]
                        v(:,1) = v(:,1)*this.ct.resolution.x;
                        v(:,2) = v(:,2)*this.ct.resolution.y;
                        v(:,3) = v(:,3)*this.ct.resolution.z;

                        % rotate surface
                        rotated_surface = v*rotMat_system_T;

                        % surface rendering
                        surface = patch('Faces',f,'Vertices',rotated_surface);
                        set(surface,'FaceColor',[0 0 1],'EdgeColor','none','FaceAlpha',.4);
                        lighting gouraud;

                    end

                    % plot projection matrix: coordinates at isocenter
                    plot3(this.rayPos(:,1),this.rayPos(:,2),this.rayPos(:,3),'k.');

                    % Plot matrix border of matrix at isocenter
                    for j = 1:stf(i).numOfRays

                        % Compute border for every bixels
                        targetPoint_vox_X_1 = stf(i).ray(j).targetPoint_bev(:,1) + this.bixelWidth;
                        targetPoint_vox_Y_1 = stf(i).ray(j).targetPoint_bev(:,2);
                        targetPoint_vox_Z_1 = stf(i).ray(j).targetPoint_bev(:,3) + this.bixelWidth;

                        targetPoint_vox_X_2 = stf(i).ray(j).targetPoint_bev(:,1) + this.bixelWidth;
                        targetPoint_vox_Y_2 = stf(i).ray(j).targetPoint_bev(:,2);
                        targetPoint_vox_Z_2 = stf(i).ray(j).targetPoint_bev(:,3) - this.bixelWidth;

                        targetPoint_vox_X_3 = stf(i).ray(j).targetPoint_bev(:,1) - this.bixelWidth;
                        targetPoint_vox_Y_3 = stf(i).ray(j).targetPoint_bev(:,2);
                        targetPoint_vox_Z_3 = stf(i).ray(j).targetPoint_bev(:,3) - this.bixelWidth;

                        targetPoint_vox_X_4 = stf(i).ray(j).targetPoint_bev(:,1) - this.bixelWidth;
                        targetPoint_vox_Y_4 = stf(i).ray(j).targetPoint_bev(:,2);
                        targetPoint_vox_Z_4 = stf(i).ray(j).targetPoint_bev(:,3) + this.bixelWidth;

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
                    plot3(this.voxTargetWorldCoords(:,1),this.voxTargetWorldCoords(:,2),this.voxTargetWorldCoords(:,3),'r.')
                    hold on;

                    % Rotated projection matrix at isocenter
                    isocenter_plane_coor = this.rayPos*rotMat_vectors_T;

                    % Plot isocenter plane
                    plot3(isocenter_plane_coor(:,1),isocenter_plane_coor(:,2),isocenter_plane_coor(:,3),'y.');

                    % Plot rotated bixels border.
                    for j = 1:stf(i).numOfRays
                        % Generate rotated projection target points.
                        targetPoint_vox_1_rotated = [stf(i).ray(j).targetPoint_bev(:,1) + this.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) + this.bixelWidth]*rotMat_vectors_T;
                        targetPoint_vox_2_rotated = [stf(i).ray(j).targetPoint_bev(:,1) + this.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) - this.bixelWidth]*rotMat_vectors_T;
                        targetPoint_vox_3_rotated = [stf(i).ray(j).targetPoint_bev(:,1) - this.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) - this.bixelWidth]*rotMat_vectors_T;
                        targetPoint_vox_4_rotated = [stf(i).ray(j).targetPoint_bev(:,1) - this.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) + this.bixelWidth]*rotMat_vectors_T;

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
                    xlabel 'X [mm]'
                    ylabel 'Y [mm]'
                    zlabel 'Z [mm]'
                    title 'lps coordinate system'
                    axis([-300 300 -300 300 -300 300]);
                    %pause(1);
                end
                % save total number of bixels
                stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);

                % Show progress
                if matRad_cfg.logLevel > 2
                    matRad_progress(i,length(this.gantryAngles));
                end
            end
        end

        function rayTargetPos = initializeRayTargetPosition(this,rayTargetPos,rotMat_vectors_T,SAD)
             % Save ray and target position in lps system.

             for j = 1:rayTargetPos.numOfRays
                  rayTargetPos.ray(j).rayPos      = rayTargetPos.ray(j).rayPos_bev*rotMat_vectors_T;
                  rayTargetPos.ray(j).targetPoint = rayTargetPos.ray(j).targetPoint_bev*rotMat_vectors_T;
             end
        end

        function stfElement = setSourceEnergyOnBeam(this,stfElement)
            % Abstract method to set source energy on beam. Implemented in Subclasses

            throw(MException('MATLAB:class:AbstractMember','This method is not implemented in the base class'));
        end


        function stfElement = finalizeBeam(this,stfElement)
            % Base method to set source energy on beam. Empty by default
            
        end
    end
end

