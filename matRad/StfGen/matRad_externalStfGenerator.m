classdef matRad_externalStfGenerator < matRad_StfGeneratorBase

    
    properties (Access = protected)
        rayPos
        ctEntryPoint
    end

    
    


    methods 
         function this = matRad_externalStfGenerator(pln)
            if nargin < 1
                pln = [];
            end
            this@matRad_StfGeneratorBase(pln);

            matRad_cfg = MatRad_Config.instance();
            addpath(fullfile(matRad_cfg.matRadRoot));

            if ~isfield(pln, 'propStf')
                matRad_cfg.dispError('no applicator information in pln struct');
            end
         end
    end

    methods (Access = protected)
        function initializePatientGeometry(this,ct, cst, visMode)
            initializePatientGeometry@matRad_StfGeneratorBase(this,ct, cst, visMode)
            matRad_cfg = MatRad_Config.instance;

            pln = this.pln; 

            if ~isfield(this.pln.propStf, 'isoCenter')
                matRad_cfg.dispWarning('No isocenter specified! Using center-of-mass of all targets!');
                this.pln.propStf.isoCenter = matRad_getIsoCenter(cst, ct);
            end

            if numel(this.pln.propStf.gantryAngles) ~= numel(this.pln.propStf.couchAngles)
                matRad_cfg.dispError('Inconsistent number of gantry and couch angles.');
            end

            if ~isnumeric(this.pln.propStf.bixelWidth) || this.pln.propStf.bixelWidth <= 0 || ~isfinite(this.pln.propStf.bixelWidth)
                matRad_cfg.dispError('Bixel width (spot distance) needs to be a real number [mm] larger than zero.');
            end
        end
        
        function stf = generateSourceGeometry(this,ct, cst, visMode)
            matRad_cfg = MatRad_Config.instance;
            pln = this.pln;

            % prepare structures necessary for particles
            SAD = this.machine.meta.SAD;

            % Define steering file like struct. Prellocating for speed.
            stf = struct;
    
            % loop over all angles
            for i = 1:length(pln.propStf.gantryAngles)

                % Correct for iso center position. Whit this correction Isocenter is
                % (0,0,0) [mm]
                coordsX = this.coordsX_vox*ct.resolution.x - pln.propStf.isoCenter(i,1);
                coordsY = this.coordsY_vox*ct.resolution.y - pln.propStf.isoCenter(i,2);
                coordsZ = this.coordsZ_vox*ct.resolution.z - pln.propStf.isoCenter(i,3);

                % Save meta information for treatment plan
                stf(i).gantryAngle   = pln.propStf.gantryAngles(i);
                stf(i).couchAngle    = pln.propStf.couchAngles(i);
                stf(i).bixelWidth    = pln.propStf.bixelWidth;
                stf(i).radiationMode = pln.radiationMode;
                stf(i).machine       = pln.machine;
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
                this.rayPos = unique(pln.propStf.bixelWidth*round([           coordsAtIsoCenterPlane(:,1) ...
                    zeros(size(coordsAtIsoCenterPlane,1),1) ...
                    coordsAtIsoCenterPlane(:,2)]/pln.propStf.bixelWidth),'rows');

                % pad ray position array if resolution of target voxel grid not sufficient
               maxCtResolution = max([ct.resolution.x ct.resolution.y ct.resolution.z]);
                if pln.propStf.bixelWidth < maxCtResolution
                    origRayPos = this.rayPos;
                    for j = -floor(maxCtResolution/pln.propStf.bixelWidth):floor(maxCtResolution/pln.propStf.bixelWidth)
                        for k = -floor(maxCtResolution/pln.propStf.bixelWidth):floor(maxCtResolution/pln.propStf.bixelWidth)
                            if abs(j)+abs(k)==0
                                continue;
                            end
                            this.rayPos = [this.rayPos; origRayPos(:,1)+j*pln.propStf.bixelWidth origRayPos(:,2) origRayPos(:,3)+k*pln.propStf.bixelWidth];
                        end
                    end
                end

                % remove spaces within rows of bixels for DAO
                if pln.propOpt.runDAO
                    % create single x,y,z vectors
                    x = this.rayPos(:,1);
                    y = this.rayPos(:,2);
                    z = this.rayPos(:,3);
                    uniZ = unique(z);
                    for j = 1:numel(uniZ)
                        x_loc = x(z == uniZ(j));
                        x_min = min(x_loc);
                        x_max = max(x_loc);
                        x = [x; (x_min:pln.propStf.bixelWidth:x_max)'];
                        y = [y; zeros((x_max-x_min)/pln.propStf.bixelWidth+1,1)];
                        z = [z; uniZ(j)*ones((x_max-x_min)/pln.propStf.bixelWidth+1,1)];
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
                rotMat_vectors_T = transpose(matRad_getRotationMatrix(pln.propStf.gantryAngles(i),pln.propStf.couchAngles(i)));


                stf(i).sourcePoint = stf(i).sourcePoint_bev*rotMat_vectors_T;

                stf(i) = this.initializeRayTargetPosition(stf(i),rotMat_vectors_T,SAD);
                
                % loop over all rays to determine meta information for each ray
                stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);

                stf(i) = this.initializeEnergy(stf(i),ct);

                if ~isfield(stf(i).ray,'energy')
                    matRad_cfg.dispError('Error generating stf struct: no suitable energies found. Check if bixelwidth is too large.');
                end
                % store total number of rays for beam-i
                stf(i).numOfRays = size(stf(i).ray,2);

                stf(i).longitudinalSpotSpacing = []; % ensures the stf(i) structs have the same structures when the next line comes
                stf(i) = this.initializePostProcessing(stf(i)); % Post Processing for Ions 

                % save total number of bixels
                stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);

                % Show progress
                if matRad_cfg.logLevel > 2
                    matRad_progress(i,length(pln.propStf.gantryAngles));
                end

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
                    plot3(this.rayPos(:,1),this.rayPos(:,2),this.rayPos(:,3),'k.');

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
                    isocenter_plane_coor = this.rayPos*rotMat_vectors_T;

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
                % save total number of bixels
                stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);

                % Show progress
                if matRad_cfg.logLevel > 2
                    matRad_progress(i,length(pln.propStf.gantryAngles));
                end
            end
        end

        function rayTargetPos = initializeRayTargetPosition(this,rayTargetPos,rotMat_vectors_T,SAD)
             % Save ray and target position in lps system.
             for j = 1:rayTargetPos.numOfRays
                  rayTargetPos.ray(j).rayPos      = rayTargetPos.ray(j).rayPos_bev*rotMat_vectors_T;
                  rayTargetPos.ray(j).targetPoint = rayTargetPos.ray(j).targetPoint_bev*rotMat_vectors_T;
                  rayTargetPos = this.initializePhotonRayPos(rayTargetPos,rotMat_vectors_T,SAD);
             end
        end

        function rays = initializeEnergy(this,rays)
        end


        function postProc = initializePostProcessing(this,postProc)
            
        end

        function photonRayPos = initializePhotonRayPos(this,photonRayPos)
        end

    end
end

