classdef (Abstract) matRad_PencilBeamEngineAbstract < DoseEngines.matRad_DoseEngineBase
    % matRad_PencilBeamEngineAbstract: abstract superclass for all dose calculation engines which are based on
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

        geometricLateralCutOff; %lateral geometric cut-off id mm, used for raytracing and geometry
        dosimetricLateralCutOff; %relative dosimetric cut-off (in fraction of values calculated)

        ssdDensityThreshold;        % Threshold for SSD computation
        useGivenEqDensityCube;      % Use the given density cube ct.cube and omit conversion from cubeHU.
        ignoreOutsideDensities;     % Ignore densities outside of cst contours
    end

    properties (SetAccess = protected)
        effectiveLateralCutOff;         %internal cutoff to be used, computed from machine/pencil-beam kernel properties and geometric/dosimetric cutoff settings
    end

    properties (SetAccess = protected, GetAccess = public)
        tmpMatrixContainers;    % temporary containers for 

        bixelsPerBeam;          % number of bixel per energy beam

        radDepthCubes;          % only stored if property set accordingly
        rotMat_system_T;        % rotation matrix for current beam
        geoDistVdoseGrid;       % geometric distance in dose grid for current beam
        rot_coordsVdoseGrid;    % Rotate coordinates for gantry movement for current beam
        radDepthVdoseGrid;      % grid for radiologica depth cube for current beam
        radDepthCube;           % radiological depth cube for current beam

        cube;                   % relative electron density / stopping power cube
        hlut;                   % hounsfield lookup table to craete relative electron density cube        
    end

    methods
        function this = matRad_PencilBeamEngineAbstract(pln)
            this = this@DoseEngines.matRad_DoseEngineBase(pln);            
        end

        function setDefaults(this)
            setDefaults@DoseEngines.matRad_DoseEngineBase(this);
            
            matRad_cfg = MatRad_Config.instance();

            %Set defaults
            this.geometricLateralCutOff       = matRad_cfg.propDoseCalc.defaultGeometricLateralCutOff;
            this.dosimetricLateralCutOff      = matRad_cfg.propDoseCalc.defaultDosimetricLateralCutOff;
            this.useGivenEqDensityCube        = matRad_cfg.propDoseCalc.defaultUseGivenEqDensityCube;
            this.ignoreOutsideDensities       = matRad_cfg.propDoseCalc.defaultIgnoreOutsideDensities;
            this.ssdDensityThreshold          = matRad_cfg.propDoseCalc.defaultSsdDensityThreshold;
        end
    end

    % Should be abstract methods but in order to satisfy the compatibility
    % with OCTAVE we can't use abstract methods. If OCTAVE at some point
    % in the far future implements this feature this should be abstract again.
    methods (Access = protected) %Abstract
        function indices = applyDoseCutOff(this)
            error('Abstract Function. Needs to be implemented!');
        end
    end

    methods (Access = protected)

        function [dij,ct,cst,stf] = calcDoseInit(this,ct,cst,stf)
            % modified inherited method of the superclass DoseEngine,
            % containing intialization which are specificly needed for
            % pencil beam calculation and not for other engines

            [dij,ct,cst,stf] = calcDoseInit@DoseEngines.matRad_DoseEngineBase(this,ct,cst,stf);
            
            % calculate rED or rSP from HU
            % Maybe we can avoid duplicating the CT here?
            if this.useGivenEqDensityCube
                matRad_cfg.dispInfo('Omitting HU to rED/rSP conversion and using existing ct.cube!\n');
            else
                ct = matRad_calcWaterEqD(ct, stf);
                %this.cube = ct.cube;
                this.hlut = ct.hlut;
            end

            %If we want to omit HU conversion check if we have a ct.cube ready
            if this.useGivenEqDensityCube && ~isfield(ct,'cube')
                matRad_cfg.dispWarning('HU Conversion requested to be omitted but no ct.cube exists! Will override and do the conversion anyway!');
                this.useGivenEqDensityCube = false;
            end

            % ignore densities outside of contours
            if this.ignoreOutsideDensities
                eraseCtDensMask = ones(prod(ct.cubeDim),1);
                eraseCtDensMask(this.VctGrid) = 0;
                for i = 1:ct.numOfCtScen
                    ct.cube{i}(eraseCtDensMask == 1) = 0;
                end
            end

            % compute SSDs
            stf = matRad_computeSSD(stf,ct,'densityThreshold',this.ssdDensityThreshold);

            % Allocate memory for quantity containers
            dij = this.allocateQuantityMatrixContainers(dij,{'physicalDose'});            
        end

        function dij = allocateQuantityMatrixContainers(this,dij,names)
            for n = 1:numel(names)
                this.tmpMatrixContainers.(names{n}) = cell(this.numOfBixelsContainer,dij.numOfScenarios);
                for i = 1:dij.numOfScenarios
                    dij.(names{n}){i} = spalloc(dij.doseGrid.numOfVoxels,this.numOfColumnsDij,1);
                end
            end

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
            if dij.numOfBeams > 1
                matRad_cfg.dispInfo('Beam %d of %d:\n',i,dij.numOfBeams);
            end

            %Reinitialize Progress:
            matRad_progress(1,1000);

            % remember beam and bixel number
            if this.calcDoseDirect
                dij.beamNum(i)    = i;
                dij.rayNum(i)     = i;
                dij.bixelNum(i)   = i;
            end

            this.bixelsPerBeam = 0;

            % convert voxel indices to real coordinates using iso center of beam i
            xCoordsV       = this.xCoordsV_vox(:)*dij.ctGrid.resolution.x-stf(i).isoCenter(1);
            yCoordsV       = this.yCoordsV_vox(:)*dij.ctGrid.resolution.y-stf(i).isoCenter(2);
            zCoordsV       = this.zCoordsV_vox(:)*dij.ctGrid.resolution.z-stf(i).isoCenter(3);
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
            %this.radDepthVdoseGrid = matRad_interpRadDepth...
            %    (ct,1,this.VctGrid,this.VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,radDepthVctGrid);

            this.radDepthVdoseGrid = this.interpRadDepth(ct,1,this.VctGrid,this.VdoseGrid,dij.ctGrid,dij.doseGrid,radDepthVctGrid);

            matRad_cfg.dispInfo('done.\n');

            %Keep rad depth cube if desired
            if this.keepRadDepthCubes
                this.radDepthCubes{i} = this.radDepthCube;
            end

            % limit rotated coordinates to positions where ray tracing is availabe
            this.rot_coordsVdoseGrid = rot_coordsVdoseGrid(~isnan(this.radDepthVdoseGrid{1}),:);

        end
        
        function radDepthVdoseGrid = interpRadDepth(~,ct,ctScenNum,V,Vcoarse,ctGrid,doseGrid,radDepthVctGrid)            
            radDepthCube                = NaN*ones(ct.cubeDim);
            radDepthCube(V(~isnan(radDepthVctGrid{1}))) = radDepthVctGrid{ctScenNum}(~isnan(radDepthVctGrid{1}));
            
            % interpolate cube - cube is now stored in Y X Z 
            coarseRadDepthCube          = matRad_interp3(ctGrid.x,ctGrid.y',ctGrid.z,radDepthCube,doseGrid.x,doseGrid.y',doseGrid.z);
            radDepthVdoseGrid{ctScenNum}  = coarseRadDepthCube(Vcoarse);
        end
        
        function ray = computeRayGeometry(this,ray,dij)
            if ~isfield(ray,'sourcePoint_bev')
                ray.sourcePoint_bev = ray.targetPoint_bev + 2*(ray.rayPos_bev - ray.targetPoint_bev);
            end

            lateralRayCutOff = this.getLateralDistanceFromDoseCutOffOnRay(ray);
            
            % Ray tracing for beam i and ray j
            [ray.ix,ray.radialDist_sq,ray.isoLatDistsX,ray.isoLatDistsZ,ray.latDistsX,ray.latDistsZ] = this.calcGeoDists(this.rot_coordsVdoseGrid, ...
                ray.sourcePoint_bev, ...
                ray.targetPoint_bev, ...
                this.machine.meta.SAD, ...
                find(~isnan(this.radDepthVdoseGrid{1})), ...
                lateralRayCutOff);

            ray.radDepths = this.radDepthVdoseGrid{1}(ray.ix);
        end
        
        function lateralRayCutOff = getLateralDistanceFromDoseCutOffOnRay(this,ray)
            lateralRayCutOff = this.effectiveLateralCutOff;
        end

        function dij = fillDij(this,dij,stf,currBeamIdx,currRayIdx,currBixelIdx,counter)
            % method for filling the dij struct with the computed dose cube
            % last step in dose calculation
            % Needs to be implemented in non abstract subclasses.
            
            if mod(counter,this.numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
                if ~this.calcDoseDirect
                    if nargin < 7
                        matRad_cfg = MatRad_Config.instance();
                        matRad_cfg.dispError('Total bixel counter not provided in fillDij');
                    end
    
                    dijColIx = (ceil(counter/this.numOfBixelsContainer)-1)*this.numOfBixelsContainer+1:counter;
                    containerIx = 1:mod(counter-1,this.numOfBixelsContainer)+1;
                    weight = 1;
                else
                    dijColIx = currBeamIdx;
                    containerIx = 1;
                    if isfield(stf(1).ray(1),'weight') && numel(stf(currBeamIdx).ray(currRayIdx).weight) >= currBixelIdx
                        weight = stf(currBeamIdx).ray(currRayIdx).weight(currBixelIdx);
                    else
                        matRad_cfg = MatRad_Config.instance();
                        matRad_cfg.dispError('No weight available for beam %d, ray %d, bixel %d',currBeamIdx,currRayIdx,currBixelIdx);
                    end
                end
                
                % Iterate through all quantities
                names = fieldnames(this.tmpMatrixContainers);
                for q = 1:numel(names)
                    qName = names{q};
                    if ~this.calcDoseDirect
                        dij.(qName){1}(:,dijColIx) = [this.tmpMatrixContainers.(qName){containerIx,1}];
                    else
                        dij.(qName){1}(:,dijColIx) = dij.(qName){1}(:,dijColIx) + weight * this.tmpMatrixContainers.(qName){containerIx,1};
                    end
                end
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

