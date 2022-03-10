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
    
    properties
       keepRadDepthCubes = false; 
    end
    
    properties (SetAccess = protected, GetAccess = public)
        
        pbCalcMode; % fine sampling mode
        
        doseTmpContainer;   % temporary container for dose calculation results
  
        effectiveLateralCutoff; % lateral cutoff for raytracing and geo calculations
        bixelsPerBeam;  % number of bixel per energy beam
                
        radDepthCubes; %only stored if property set accordingly)
    end
    
    properties (Access = protected)
        rotMat_system_T; % rotation matrix for current beam
        geoDistVdoseGrid;   % geometric distance in dose grid for current beam
        rot_coordsVdoseGrid;    % Rotate coordinates for gantry movement for current beam
        radDepthVdoseGrid;  % grid for radiologica depth cube for current beam
        radDepthCube;   % radiological depth cube for current beam
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
            % Method for initializing the beams for analytical pencil beam
            % dose calculation
            %
            % call
            %   obj.calcDoseInitBeam(ct,stf,dij,i)
            %
            % input
            %   ct:                         matRad ct struct
            %   stf:                        matRad steering information struct
            %   dij:                        matRad dij struct
            %   i:                          index of beam
            %
            % output
            %   dij:                        updated dij struct
            
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

            obj.rotMat_system_T = matRad_getRotationMatrix(stf(i).gantryAngle,stf(i).couchAngle);

            % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
            rot_coordsV         = coordsV*obj.rotMat_system_T;
            rot_coordsVdoseGrid = coordsVdoseGrid*obj.rotMat_system_T;

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
            if strcmp(obj.pbCalcMode, 'fineSampling') || obj.keepRadDepthCubes
                [radDepthVctGrid, obj.radDepthCube] = matRad_rayTracing(stf(i),ct,obj.VctGrid,rot_coordsV,obj.effectiveLateralCutoff);
                obj.radDepthCube{1} = matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y,   dij.ctGrid.z, obj.radDepthCube{1}, ...
                                                dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest');
            else
                radDepthVctGrid = matRad_rayTracing(stf(i),ct,obj.VctGrid,rot_coordsV,obj.effectiveLateralCutoff);
            end
            
            % interpolate radiological depth cube to dose grid resolution
            obj.radDepthVdoseGrid = matRad_interpRadDepth...
                (ct,1,obj.VctGrid,obj.VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,radDepthVctGrid);
            
            matRad_cfg.dispInfo('done.\n');

            %Keep rad depth cube if desired
            if obj.keepRadDepthCubes
                obj.radDepthCubes{i} = obj.radDepthCube;
            end

            % limit rotated coordinates to positions where ray tracing is availabe
            obj.rot_coordsVdoseGrid = rot_coordsVdoseGrid(~isnan(obj.radDepthVdoseGrid{1}),:);
            
        end
        
    end
    
    methods (Static)
        
        function [ix,rad_distancesSq,isoLatDistsX,isoLatDistsZ,latDistsX,latDistsZ] = ...
              calcGeoDists(rot_coords_bev, sourcePoint_bev, targetPoint_bev, SAD, radDepthIx, lateralCutOff)
            % matRad calculation of lateral distances from central ray 
            % used for dose calculation
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
        
        function [ixNew,bixelDoseNew] =  dijSampling(ix,bixelDose,radDepthV,rad_distancesSq,sType,Param)
            % matRad dij sampling function 
            % This function samples. 
            % 
            % call
            %   [ixNew,bixelDoseNew] =
            %   obj.dijSampling(ix,bixelDose,radDepthV,rad_distancesSq,sType,Param)
            %
            % input
            %   ix:               indices of voxels where we want to compute dose influence data
            %   bixelDose:        dose at specified locations as linear vector
            %   radDepthV:        radiological depth vector
            %   rad_distancesSq:  squared radial distance to the central ray
            %   sType:            can either be set to 'radius' or 'dose'. These are two different ways 
            %                     to determine dose values that are keept as they are and dose values used for sampling
            %   Param:            In the case of radius based sampling, dose values having a radial 
            %                     distance below r0 [mm] are keept anyway and sampling is only done beyond r0. 
            %                     In the case of dose based sampling, dose values having a relative dose greater 
            %                     the threshold [0...1] are keept and sampling is done for dose values below the relative threshold  
            %
            % output
            %   ixNew:            reduced indices of voxels where we want to compute dose influence data
            %   bixelDoseNew      reduced dose at specified locations as linear vector
            %
            % References
            %   [1] http://dx.doi.org/10.1118/1.1469633
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Copyright 2016 the matRad development team. 
            % 
            % This file is part of the matRad project. It is subject to the license 
            % terms in the LICENSE file found in the top-level directory of this 
            % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
            % of the matRad project, including this file, may be copied, modified, 
            % propagated, or distributed except according to the terms contained in the 
            % LICENSE file.
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% define default parameters as a fallback 
            defaultType                = 'radius';
            deltaRadDepth              = 5;                       % step size of radiological depth
            defaultLatCutOff           = 25;                      % default lateral cut off
            defaultrelDoseThreshold    = 0.01;                    % default relative dose threshold

            relDoseThreshold           = defaultrelDoseThreshold;
            LatCutOff                  = defaultLatCutOff;
            Type                       = sType;

            % if the input index vector is of type logical convert it to linear indices
            if islogical(ix)
               ix = find(ix); 
            end

            %% parse inputs
            if sum(strcmp(sType,{'radius','dose'})) == 0
               Type = defaultType;
            end

            % if an parameter is provided then use it
            if nargin>5   
                if exist('Param','var')
                     if strcmp(sType,'radius')
                       LatCutOff = Param;
                    elseif strcmp(sType,'dose')
                       relDoseThreshold = Param;
                    end
                end
            end

            %% remember dose values inside the inner core
            switch  Type
                case {'radius'}
                ixCore      = rad_distancesSq < LatCutOff^2;                 % get voxels indices having a smaller radial distance than r0
                case {'dose'}
                ixCore      = bixelDose > relDoseThreshold * max(bixelDose); % get voxels indices having a greater dose than the thresholdDose
            end

            bixelDoseCore       = bixelDose(ixCore);                         % save dose values that are not affected by sampling

            if all(ixCore)
                %% all bixels are in the core
                %exit function with core dose only
                ixNew = ix;
                bixelDoseNew = bixelDoseCore;
            else
                logIxTail           = ~ixCore;                                   % get voxels indices beyond r0
                linIxTail           = find(logIxTail);                           % convert logical index to linear index
                numTail             = numel(linIxTail);
                bixelDoseTail       = bixelDose(linIxTail);                      % dose values that are going to be reduced by sampling
                ixTail              = ix(linIxTail);                             % indices that are going to be reduced by sampling

                %% sample for each radiological depth the lateral halo dose
                radDepthTail        = (radDepthV(linIxTail));                    % get radiological depth in the tail

                % cluster radiological dephts to reduce computations
                B_r                 = int32(ceil(radDepthTail));                 % cluster radiological depths;
                maxRadDepth         = double(max(B_r));
                C                   = int32(linspace(0,maxRadDepth,round(maxRadDepth)/deltaRadDepth));     % coarse clustering of rad depths

                ixNew               = zeros(numTail,1);                          % inizialize new index vector
                bixelDoseNew        = zeros(numTail,1);                          % inizialize new dose vector
                linIx               = int32(1:1:numTail)';
                IxCnt               = 1;

                %% loop over clustered radiological depths
                for i = 1:numel(C)-1
                    ixTmp              = linIx(B_r >= C(i) & B_r < C(i+1));      % extracting sub indices
                    if isempty(ixTmp)
                        continue
                    end
                    subDose            = bixelDoseTail(ixTmp);                   % get tail dose in current cluster
                    subIx              = ixTail(ixTmp);                          % get indices in current cluster
                    thresholdDose      = max(subDose);
                    r                  = rand(numel(subDose),1);                 % get random samples
                    ixSamp             = r<=(subDose/thresholdDose);
                    NumSamples         = sum(ixSamp);

                    ixNew(IxCnt:IxCnt+NumSamples-1,1)        = subIx(ixSamp);    % save new indices
                    bixelDoseNew(IxCnt:IxCnt+NumSamples-1,1) = thresholdDose;    % set the dose
                    IxCnt = IxCnt + NumSamples;
                end


                % cut new vectors and add inner core values
                ixNew        = [ix(ixCore);    ixNew(1:IxCnt-1)];
                bixelDoseNew = [bixelDoseCore; bixelDoseNew(1:IxCnt-1)];
            end
            
        end
        
    end
   
end

