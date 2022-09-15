classdef (Abstract) matRad_DoseEnginePencilBeam < DoseEngines.matRad_DoseEngine
    % matRad_DoseEnginePencilBeam: abstract superclass for all dose calculation engines which are based on
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

        enableDijSampling = true;
        dijSampling;             %struct with lateral dij sampling parameters

        geometricLateralCutOff; %lateral geometric cut-off id mm, used for raytracing and geometry
        dosimetricLateralCutOff; %relative dosimetric cut-off (in fraction of values calculated)
    end

    properties (SetAccess = protected)
        effectiveLateralCutOff;         %internal cutoff to be used, computed from machine/pencil-beam kernel properties and geometric/dosimetric cutoff settings
    end

    properties (SetAccess = protected, GetAccess = public)
        doseTmpContainer;       % temporary container for dose calculation results

        bixelsPerBeam;          % number of bixel per energy beam

        radDepthCubes;          % only stored if property set accordingly
        rotMat_system_T;        % rotation matrix for current beam
        geoDistVdoseGrid;       % geometric distance in dose grid for current beam
        rot_coordsVdoseGrid;    % Rotate coordinates for gantry movement for current beam
        radDepthVdoseGrid;      % grid for radiologica depth cube for current beam
        radDepthCube;           % radiological depth cube for current beam
    end

    methods
        function this = matRad_DoseEnginePencilBeam()
            this = this@DoseEngines.matRad_DoseEngine();
            
            matRad_cfg = MatRad_Config.instance();

            %Set defaults
            this.geometricLateralCutOff       = matRad_cfg.propDoseCalc.defaultGeometricLateralCutOff;
            this.dosimetricLateralCutOff      = matRad_cfg.propDoseCalc.defaultDosimetricLateralCutOff;

            %Set additional parameters
            this.dijSampling.relDoseThreshold = 0.01;
            this.dijSampling.latCutOff        = 20;
            this.dijSampling.type             = 'radius';
            this.dijSampling.deltaRadDepth    = 5;
        end
    end

    % Should be abstract methods but in order to satisfy the compatibility
    % with OCTAVE we can't use abstract methods. If OCTAVE at some point
    % in the far future implements this feature this should be abstract again.
    methods (Access = protected) %Abstract


        function dij = fillDij(this,dij,stf,counter)
            % method for filling the dij struct with the computed dose cube
            % last step in dose calculation
            % Needs to be implemented in non abstract subclasses.
            error('Funktion needs to be implemented!');
        end


        function dij = fillDijDirect(this,dij,stf,currBeamIdx,currRayIdx,currBixelIdx)
            % method for filling the dij struct, when using a direct dose
            % calcultion
            % Needs to be implemented in non abstract subclasses,
            % when direct calc shoulb be utilizable.
            error('Funktion needs to be implemented!');
        end

    end

    methods (Access = protected)

        function [dij,ct,cst,stf,pln] = calcDoseInit(this,ct,cst,stf,pln)
            % modified inherited method of the superclass DoseEngine,
            % containing intialization which are specificly needed for
            % pencil beam calculation and not for other engines

            [dij,ct,cst,stf,pln] = calcDoseInit@DoseEngines.matRad_DoseEngine(this,ct,cst,stf,pln);

            % Allocate memory for dose_temp cell array
            this.doseTmpContainer     = cell(this.numOfBixelsContainer,dij.numOfScenarios);
        end

        function dij = calcDoseInitBeam(this,dij,ct,cst,stf,i)
            % Method for initializing the beams for analytical pencil beam
            % dose calculation
            %
            % call
            %   this.calcDoseInitBeam(ct,stf,dij,i)
            %
            % input
            %   ct:                         matRad ct struct
            %   cst:                        matRad cst struct
            %   stf:                        matRad steering information struct
            %   i:                          index of beam
            %
            % output
            %   dij:                        updated dij struct

            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispInfo('Beam %d of %d:\n',i,dij.numOfBeams);

            % remember beam and bixel number
            if this.calcDoseDirect
                dij.beamNum(i)    = i;
                dij.rayNum(i)     = i;
                dij.bixelNum(i)   = i;
            end

            this.bixelsPerBeam = 0;

            % convert voxel indices to real coordinates using iso center of beam i
            xCoordsV       = this.xCoordsV_vox(:)*ct.resolution.x-stf(i).isoCenter(1);
            yCoordsV       = this.yCoordsV_vox(:)*ct.resolution.y-stf(i).isoCenter(2);
            zCoordsV       = this.zCoordsV_vox(:)*ct.resolution.z-stf(i).isoCenter(3);
            coordsV        = [xCoordsV yCoordsV zCoordsV];

            xCoordsVdoseGrid = this.xCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.x-stf(i).isoCenter(1);
            yCoordsVdoseGrid = this.yCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.y-stf(i).isoCenter(2);
            zCoordsVdoseGrid = this.zCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.z-stf(i).isoCenter(3);
            coordsVdoseGrid  = [xCoordsVdoseGrid yCoordsVdoseGrid zCoordsVdoseGrid];

            % Get Rotation Matrix
            % Do not transpose matrix since we usage of row vectors &
            % transformation of the coordinate system need double transpose

            this.rotMat_system_T = matRad_getRotationMatrix(stf(i).gantryAngle,stf(i).couchAngle);

            % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
            rot_coordsV         = coordsV*this.rotMat_system_T;
            rot_coordsVdoseGrid = coordsVdoseGrid*this.rotMat_system_T;

            rot_coordsV(:,1) = rot_coordsV(:,1)-stf(i).sourcePoint_bev(1);
            rot_coordsV(:,2) = rot_coordsV(:,2)-stf(i).sourcePoint_bev(2);
            rot_coordsV(:,3) = rot_coordsV(:,3)-stf(i).sourcePoint_bev(3);

            rot_coordsVdoseGrid(:,1) = rot_coordsVdoseGrid(:,1)-stf(i).sourcePoint_bev(1);
            rot_coordsVdoseGrid(:,2) = rot_coordsVdoseGrid(:,2)-stf(i).sourcePoint_bev(2);
            rot_coordsVdoseGrid(:,3) = rot_coordsVdoseGrid(:,3)-stf(i).sourcePoint_bev(3);

            % calculate geometric distances
            this.geoDistVdoseGrid{1}= sqrt(sum(rot_coordsVdoseGrid.^2,2));
            % Calculate radiological depth cube
            matRad_cfg.dispInfo('matRad: calculate radiological depth cube... ');
            if this.keepRadDepthCubes
                [radDepthVctGrid, this.radDepthCube] = matRad_rayTracing(stf(i),ct,this.VctGrid,rot_coordsV,this.effectiveLateralCutOff);
                this.radDepthCube{1} = matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y,   dij.ctGrid.z, this.radDepthCube{1}, ...
                    dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest');
            else
                radDepthVctGrid = matRad_rayTracing(stf(i),ct,this.VctGrid,rot_coordsV,this.effectiveLateralCutOff);
            end

            % interpolate radiological depth cube to dose grid resolution
            this.radDepthVdoseGrid = matRad_interpRadDepth...
                (ct,1,this.VctGrid,this.VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,radDepthVctGrid);

            matRad_cfg.dispInfo('done.\n');

            %Keep rad depth cube if desired
            if this.keepRadDepthCubes
                this.radDepthCubes{i} = this.radDepthCube;
            end

            % limit rotated coordinates to positions where ray tracing is availabe
            this.rot_coordsVdoseGrid = rot_coordsVdoseGrid(~isnan(this.radDepthVdoseGrid{1}),:);

        end

        function [ixNew,bixelDoseNew] =  sampleDij(this,ix,bixelDose,radDepthV,rad_distancesSq,bixelWidth)
            % matRad dij sampling function
            % This function samples.
            %
            % call
            %   [ixNew,bixelDoseNew] =
            %   this.sampleDij(ix,bixelDose,radDepthV,rad_distancesSq,sType,Param)
            %
            % input
            %   ix:               indices of voxels where we want to compute dose influence data
            %   bixelDose:        dose at specified locations as linear vector
            %   radDepthV:        radiological depth vector
            %   rad_distancesSq:  squared radial distance to the central ray
            %   bixelWidth:       bixelWidth as set in pln (optional)
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

            relDoseThreshold           = this.dijSampling.relDoseThreshold;
            LatCutOff                  = this.dijSampling.latCutOff;
            Type                       = this.dijSampling.type;
            deltaRadDepth              = this.dijSampling.deltaRadDepth;

            % if the input index vector is of type logical convert it to linear indices
            if islogical(ix)
                ix = find(ix);
            end

            %Increase sample cut-off by bixel width if given
            if nargin == 6 && ~isempty(bixelWidth)
                LatCutOff = LatCutOff + bixelWidth;
            end

            %% remember dose values inside the inner core
            switch  Type
                case 'radius'
                    ixCore      = rad_distancesSq < LatCutOff^2;                 % get voxels indices having a smaller radial distance than r0
                case 'dose'
                    ixCore      = bixelDose > relDoseThreshold * max(bixelDose); % get voxels indices having a greater dose than the thresholdDose
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Dij Sampling mode ''%s'' not known!',Type);
            end

            bixelDoseCore = bixelDose(ixCore);                         % save dose values that are not affected by sampling

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

    methods (Static)

        function [ix,rad_distancesSq,isoLatDistsX,isoLatDistsZ,latDistsX,latDistsZ] = ...
                calcGeoDists(rot_coords_bev, sourcePoint_bev, targetPoint_bev, SAD, radDepthIx, lateralCutOff)
            % matRad calculation of lateral distances from central ray
            % used for dose calculation
            %
            % call
            %   [ix,rad_distancesSq,isoLatDistsX,isoLatDistsZ] = ...
            %           this.calcGeoDists(rot_coords_bev, ...
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

    %% deprecated properties
    properties (Dependent)
        geometricCutOff;                %deprecated property, replaced with geometricLateralCutOff
    end

    methods
        function set.geometricCutOff(this,geoCutOff)
            this.geometricLateralCutOff = geoCutOff;
            this.warnDeprecatedEngineProperty('geometricCutOff','','geometricLateralCutOff');
        end
        function geoCutOff = get.geometricCutOff(this)
            geoCutOff = this.geometricLateralCutOff;
            this.warnDeprecatedEngineProperty('geometricCutOff','','geometricLateralCutOff');
        end
    end

end

