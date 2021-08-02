classdef (Abstract) matRad_AnalyticalPencilBeamEngine < DoseEngines.matRad_DoseEngine
    % matRad_AnalyticalPencilBeamEngine: abstract superclass for all dose calculation engines which are based on 
    %   analytical pencil beam calculation 
    %   for more informations see superclass
    %   DoseEngines.matRad_DoseEngine
    %   MatRad_Config MatRad Configuration class
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2019 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (SetAccess = protected, GetAccess = public)
        
        pbCalcMode; % fine sampling mode
        
        doseTmpContainer;   % temporary container for dose calculation results
          
        geoDistVdoseGrid;   % geometric distance in dose grid
        rot_coordsVdoseGrid;    % Rotate coordinates for gantry movement
        radDepthsMat;   % radiological depth cube container
        radDepthVdoseGrid;  % grid for radiologica depth cube
       
        effectiveLateralCutoff; % lateral cutoff for raytracing and geo calculations
        
        bixelsPerBeam;  % number of bixel per energy beam
        
    end
    
    % Should be abstract methods but in order to satisfy the compatibility
    % with OCTAVE we can't use abstract methods. If OCTAVE at some point 
    % in the far future implements this feature this should be abstract again.
    methods (Access = protected) %Abstract

        
        function dij = fillDij(obj,dij,stf,counter) 
        % method for filling the dij struct with the computed dose cube
        % last step in dose calculation
        % Needs to be implemented in non abstract subclasses.  
            error('Funktion needs to be implemented!');
        end

        
        function dij = fillDijDirect(obj,dij,stf,currBeamIdx,currRayIdx,currBixelIdx)
        % method for filling the dij struct, when using a direct dose
        % calcultion
        % Needs to be implemented in non abstract subclasses, 
        % when direct calc shoulb be utilizable.    
            error('Funktion needs to be implemented!');
        end

    end
    
    methods (Access = protected)
        
        function [ct,stf,pln,dij] = calcDoseInit(obj,ct,stf,pln,cst)
            % modified inherited method of the superclass DoseEngine,
            % containing intialization which are specificly needed for
            % pencil beam calculation and not for other engines
            
            [ct,stf,pln,dij] = calcDoseInit@DoseEngines.matRad_DoseEngine(obj,ct,stf,pln,cst);
            
            % Allocate memory for dose_temp cell array
            obj.doseTmpContainer     = cell(obj.numOfBixelsContainer,dij.numOfScenarios);
            
            
            
        end
        
        function dij = calcDoseInitBeam(obj,ct,stf,dij,i)
            
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispInfo('Beam %d of %d:\n',i,dij.numOfBeams);

            % remember beam and bixel number
            if obj.calcDoseDirect
                dij.beamNum(i)    = i;
                dij.rayNum(i)     = i;
                dij.bixelNum(i)   = i;
            end

            obj.bixelsPerBeam = 0;

            % convert voxel indices to real coordinates using iso center of beam i
            xCoordsV       = obj.xCoordsV_vox(:)*ct.resolution.x-stf(i).isoCenter(1);
            yCoordsV       = obj.yCoordsV_vox(:)*ct.resolution.y-stf(i).isoCenter(2);
            zCoordsV       = obj.zCoordsV_vox(:)*ct.resolution.z-stf(i).isoCenter(3);
            coordsV        = [xCoordsV yCoordsV zCoordsV];

            xCoordsVdoseGrid = obj.xCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.x-stf(i).isoCenter(1);
            yCoordsVdoseGrid = obj.yCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.y-stf(i).isoCenter(2);
            zCoordsVdoseGrid = obj.zCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.z-stf(i).isoCenter(3);
            coordsVdoseGrid  = [xCoordsVdoseGrid yCoordsVdoseGrid zCoordsVdoseGrid];

            % Get Rotation Matrix
            % Do not transpose matrix since we usage of row vectors &
            % transformation of the coordinate system need double transpose

            rotMat_system_T = matRad_getRotationMatrix(stf(i).gantryAngle,stf(i).couchAngle);

            % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
            rot_coordsV         = coordsV*rotMat_system_T;
            rot_coordsVdoseGrid = coordsVdoseGrid*rotMat_system_T;

            rot_coordsV(:,1) = rot_coordsV(:,1)-stf(i).sourcePoint_bev(1);
            rot_coordsV(:,2) = rot_coordsV(:,2)-stf(i).sourcePoint_bev(2);
            rot_coordsV(:,3) = rot_coordsV(:,3)-stf(i).sourcePoint_bev(3);

            rot_coordsVdoseGrid(:,1) = rot_coordsVdoseGrid(:,1)-stf(i).sourcePoint_bev(1);
            rot_coordsVdoseGrid(:,2) = rot_coordsVdoseGrid(:,2)-stf(i).sourcePoint_bev(2);
            rot_coordsVdoseGrid(:,3) = rot_coordsVdoseGrid(:,3)-stf(i).sourcePoint_bev(3);

            % calculate geometric distances
            obj.geoDistVdoseGrid{1}= sqrt(sum(rot_coordsVdoseGrid.^2,2));
            % Calculate radiological depth cube
            matRad_cfg.dispInfo('matRad: calculate radiological depth cube... ');
            if strcmp(obj.pbCalcMode, 'fineSampling')
                [radDepthVctGrid, obj.radDepthsMat] = matRad_rayTracing(stf(i),ct,obj.VctGrid,rot_coordsV,obj.effectiveLateralCutoff);
            else
                radDepthVctGrid = matRad_rayTracing(stf(i),ct,obj.VctGrid,rot_coordsV,obj.effectiveLateralCutoff);
            end
            matRad_cfg.dispInfo('done.\n');

            % interpolate radiological depth cube to dose grid resolution
            obj.radDepthVdoseGrid = matRad_interpRadDepth...
                (ct,1,obj.VctGrid,obj.VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,radDepthVctGrid);

            if exist('radDepthsMat', 'var')
                % interpolate radiological depth cube used for fine sampling to dose grid resolution
                obj.radDepthsMat{1} = matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y,   dij.ctGrid.z, obj.radDepthsMat{1}, ...
                                                dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest');
            end

            % limit rotated coordinates to positions where ray tracing is availabe
            obj.rot_coordsVdoseGrid = rot_coordsVdoseGrid(~isnan(obj.radDepthVdoseGrid{1}),:);
            
        end
    end
    
    methods (Static)
        
        function [ix,rad_distancesSq,isoLatDistsX,isoLatDistsZ,latDistsX,latDistsZ] = ...
              calcGeoDists(rot_coords_bev, sourcePoint_bev, targetPoint_bev, SAD, radDepthIx, lateralCutOff)
            % matRad calculation of lateral distances from central ray 
            % used for dose calcultion
            % 
            % call
            %   [ix,rad_distancesSq,isoLatDistsX,isoLatDistsZ] = ...
            %           obj.calcGeoDists(rot_coords_bev, ...
            %                               sourcePoint_bev, ...
            %                               targetPoint_bev, ...
            %                               SAD, ...
            %                               radDepthIx, ...
            %                               lateralCutOff)
            %
            % input
            %   rot_coords_bev:     coordinates in bev of the voxels with index V,
            %                       where also ray tracing results are availabe 
            %   sourcePoint_bev:    source point in voxel coordinates in beam's eye view
            %   targetPoint_bev:    target point in voxel coordinated in beam's eye view
            %   SAD:                source-to-axis distance
            %   radDepthIx:         sub set of voxels for which radiological depth
            %                       calculations are available
            %   lateralCutOff:      lateral cutoff specifying the neighbourhood for
            %                       which dose calculations will actually be performed
            %
            % output
            %   ix:                 indices of voxels where we want to compute dose
            %                       influence data
            %   rad_distancesSq:    squared radial distance to the central ray (where the
            %                       actual computation of the radiological depth takes place)
            %   isoLatDistsX:       lateral x-distance to the central ray projected to
            %                       iso center plane
            %   isoLatDistsZ:       lateral z-distance to the central ray projected to
            %                       iso center plane
            %   latDistsX:          lateral x-distance to the central ray
            %   latDistsZ:          lateral z-distance to the central ray
            %
            %
            % References
            %   -
            %
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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ROTATE A SINGLE BEAMLET AND ALIGN WITH BEAMLET WHO PASSES THROUGH
            % ISOCENTER

            % Put [0 0 0] position in the source point for beamlet who passes through
            % isocenter
            a = -sourcePoint_bev';

            % Normalize the vector
            a = a/norm(a);

            % Put [0 0 0] position in the source point for a single beamlet
            b = (targetPoint_bev - sourcePoint_bev)';

            % Normalize the vector
            b = b/norm(b);

            % Define function for obtain rotation matrix.
            if all(a==b) % rotation matrix corresponds to eye matrix if the vectors are the same
                rot_coords_temp = rot_coords_bev;
            else
                % Define rotation matrix
                ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
                R   = eye(3) + ssc(cross(a,b)) + ssc(cross(a,b))^2*(1-dot(a,b))/(norm(cross(a,b))^2);

                % Rotate every CT voxel 
                rot_coords_temp = rot_coords_bev*R;
            end

            % Put [0 0 0] position CT in center of the beamlet.
            latDistsX = rot_coords_temp(:,1) + sourcePoint_bev(1);
            latDistsZ = rot_coords_temp(:,3) + sourcePoint_bev(3);

            % check of radial distance exceeds lateral cutoff (projected to iso center)
            rad_distancesSq = latDistsX.^2 + latDistsZ.^2;
            subsetMask = rad_distancesSq ./ rot_coords_temp(:,2).^2 <= lateralCutOff^2 /SAD^2;

            % return index list within considered voxels
            ix = radDepthIx(subsetMask);

            % return radial distances squared
            rad_distancesSq = rad_distancesSq(subsetMask);

            % return x & z distance
            % if nargout > 2
            %    isoLatDistsX = latDistsX(subsetMask)./rot_coords_temp(subsetMask,2)*SAD;
            %    isoLatDistsZ = latDistsZ(subsetMask)./rot_coords_temp(subsetMask,2)*SAD; 
            % end


            % latDists
            if nargout > 4
                % latDists
                latDistsX = latDistsX(subsetMask);
                latDistsZ = latDistsZ(subsetMask);
                isoLatDistsX = latDistsX./rot_coords_temp(subsetMask,2)*SAD;
                isoLatDistsZ = latDistsZ./rot_coords_temp(subsetMask,2)*SAD; 
            else
                % lateral distances projected to iso center plane
                isoLatDistsX = latDistsX(subsetMask)./rot_coords_temp(subsetMask,2)*SAD;
                isoLatDistsZ = latDistsZ(subsetMask)./rot_coords_temp(subsetMask,2)*SAD; 
            end


        end
    end
   
end

