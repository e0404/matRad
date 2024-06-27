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

        numOfDijFillSteps = 10;     % Number of times during dose calculation the temporary containers are moved to a sparse matrix
    end

    properties (SetAccess = protected)
        effectiveLateralCutOff;         %internal cutoff to be used, computed from machine/pencil-beam kernel properties and geometric/dosimetric cutoff settings
    end

    properties (SetAccess = protected, GetAccess = public)
        tmpMatrixContainers;    % temporary containers for 
        numOfBixelsContainer;   % number of used bixel container

        radDepthCubes = {};     % only stored if property set accordingly

        cubeWED;                % relative electron density / stopping power cube
        hlut;                   % hounsfield lookup table to craete relative electron density cube    
    end

    methods
        function this = matRad_PencilBeamEngineAbstract(pln)
            if nargin < 1
                pln = [];
            end

            this = this@DoseEngines.matRad_DoseEngineBase(pln);            
        end

        function setDefaults(this)
            setDefaults@DoseEngines.matRad_DoseEngineBase(this);
            
            matRad_cfg = MatRad_Config.instance();

            %Set defaults
            this.geometricLateralCutOff       = matRad_cfg.defaults.propDoseCalc.geometricLateralCutOff;
            this.dosimetricLateralCutOff      = matRad_cfg.defaults.propDoseCalc.dosimetricLateralCutOff;
            this.useGivenEqDensityCube        = matRad_cfg.defaults.propDoseCalc.useGivenEqDensityCube;
            this.ignoreOutsideDensities       = matRad_cfg.defaults.propDoseCalc.ignoreOutsideDensities;
            this.ssdDensityThreshold          = matRad_cfg.defaults.propDoseCalc.ssdDensityThreshold;
        end    
    end

    % Should be abstract methods but in order to satisfy the compatibility
    % with OCTAVE we can't use abstract methods. If OCTAVE at some point
    % in the far future implements this feature this should be abstract again.
    methods (Access = protected) %Abstract
        function bixel = computeBixel(this,currRay,k)
            throw(MException('MATLAB:class:AbstractMember','Abstract function computeBixel of your PencilBeam DoseEngine needs to be implemented!'));
        end
    end

    methods (Access = protected)
        
        function dij = calcDose(this,ct,cst,stf)
            matRad_cfg = MatRad_Config.instance();

            % initialize
            dij = this.initDoseCalc(ct,cst,stf);

            for shiftScen = 1:this.multScen.totNumShiftScen

                %Find first instance of the shift to select the shift values
                ixShiftScen = find(this.multScen.linearMask(:,2) == shiftScen,1);

                scenStf = stf;
                % manipulate isocenter
                for k = 1:numel(scenStf)
                    scenStf(k).isoCenter = scenStf(k).isoCenter + this.multScen.isoShift(ixShiftScen,:);
                end

                if this.multScen.totNumShiftScen > 1
                    matRad_cfg.dispInfo('\tShift scenario %d of %d: \n',shiftScen,this.multScen.totNumShiftScen);
                end

                bixelCounter = 0;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for i = 1:dij.numOfBeams % loop over all beams

                    %Initialize Beam Geometry
                    currBeam = this.initBeam(dij,ct,cst,scenStf,i);

                    %Keep tabs on bixels computed in this beam
                    bixelBeamCounter = 0;

                    %Ray calculation
                    for j = 1:currBeam.numOfRays % loop over all rays / for photons we only have one bixel per ray! For field based dose calc, a ray equals a shape

                        %Initialize Ray Geometry
                        currRay = this.initRay(currBeam,j);

                        for ctScen = 1:this.multScen.numOfCtScen
                            for rangeShiftScen = 1:this.multScen.totNumRangeScen
                                fullScenIdx = this.multScen.sub2scenIx(ctScen,shiftScen,rangeShiftScen);
 
                                if this.multScen.scenMask(fullScenIdx)
                                    %TODO: This shows we probably need
                                    %better scenario management
                                    %Gets linear index in scenario cell array
                                    

                                    scenRay = this.extractSingleScenarioRay(currRay,fullScenIdx);    
                                    
                                    for k = 1:currRay.numOfBixels
                                        %Bixel Computation
                                        currBixel = this.computeBixel(scenRay,k);

                                        % save computation time and memory
                                        % by sequentially filling the sparse matrix dose.dij from the cell array
                                        dij = this.fillDij(currBixel,dij,scenStf,fullScenIdx,i,j,k,bixelCounter + k);
                                    end
                                end
                            end
                        end

                        % Progress Update & Bookkeeping
                        bixelCounter = bixelCounter + currRay.numOfBixels;
                        bixelBeamCounter = bixelBeamCounter + currRay.numOfBixels;
                        this.progressUpdate(bixelCounter,dij.totalNumOfBixels);
                    end
                end
            end

            %Finalize dose calculation
            dij = this.finalizeDose(dij);
        end

        function dij = initDoseCalc(this,ct,cst,stf)
            % modified inherited method of the superclass DoseEngine,
            % containing intialization which are specificly needed for
            % pencil beam calculation and not for other engines

            dij = initDoseCalc@DoseEngines.matRad_DoseEngineBase(this,ct,cst,stf);
            
            matRad_cfg = MatRad_Config.instance();
            
            % calculate rED or rSP from HU or take provided wedCube
            if this.useGivenEqDensityCube && ~isfield(ct,'cube')
                matRad_cfg.dispWarning('HU Conversion requested to be omitted but no ct.cube exists! Will override and do the conversion anyway!');
                this.useGivenEqDensityCube = false;
            end

            if this.useGivenEqDensityCube
                matRad_cfg.dispInfo('Omitting HU to rED/rSP conversion and using existing ct.cube!\n');
            else
                ct = matRad_calcWaterEqD(ct, stf); % Maybe we can avoid duplicating the CT here?
            end

            this.cubeWED = ct.cube;
            if isfield(ct,'hlut')
                this.hlut = ct.hlut;
            end

            % ignore densities outside of contours
            if this.ignoreOutsideDensities
                eraseCtDensMask = ones(prod(ct.cubeDim),1);
                eraseCtDensMask(this.VctGrid) = 0;
                for i = 1:ct.numOfCtScen
                    this.cubeWED{i}(eraseCtDensMask == 1) = 0;
                end
            end

            % Allocate memory for quantity containers
            dij = this.allocateQuantityMatrixContainers(dij,{'physicalDose'});            
        end

        function dij = allocateQuantityMatrixContainers(this,dij,names)
            if this.calcDoseDirect
                this.numOfBixelsContainer = 1;
            else
                this.numOfBixelsContainer = ceil(dij.totalNumOfBixels/this.numOfDijFillSteps);
            end
            
            %Loop over all requested quantities
            for n = 1:numel(names)
                %Create Cell arrays for container and dij
                szContainer = [this.numOfBixelsContainer size(this.multScen.scenMask)];
                this.tmpMatrixContainers.(names{n}) = cell(szContainer);
                dij.(names{n}) = cell(size(this.multScen.scenMask));
                
                %Now preallocate a matrix in each active scenario using the
                %scenmask
                if this.calcDoseDirect
                    dij.(names{n})(this.multScen.scenMask) = {zeros(dij.doseGrid.numOfVoxels,this.numOfColumnsDij)};
                else
                    %We preallocate a sparse matrix with sparsity of
                    %1e-3 to make the filling slightly faster
                    %TODO: the preallocation could probably
                    %have more accurate estimates
                    dij.(names{n})(this.multScen.scenMask) = {spalloc(dij.doseGrid.numOfVoxels,this.numOfColumnsDij,round(prod(dij.doseGrid.numOfVoxels,this.numOfColumnsDij)*1e-3))};
                end
            end
        end

        function currBeam = initBeam(this,dij,ct,cst,stf,i)
            % Method for initializing the beams for analytical pencil beam
            % dose calculation
            %
            % call
            %   this.initBeam(ct,stf,dij,i)
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
            if numel(stf) > 1
                matRad_cfg.dispInfo('Beam %d of %d:\n',i,numel(stf));
            end

            currBeam = stf(i);
            currBeam.beamIndex = i;

            %Adjust for isocenter offset
            %change the stf
            currBeam.isoCenter = currBeam.isoCenter + dij.doseGrid.isoCenterOffset;

            % convert voxel indices to real coordinates using iso center of beam i
            xCoordsV       = this.xCoordsV_vox(:)*dij.ctGrid.resolution.x-currBeam.isoCenter(1);
            yCoordsV       = this.yCoordsV_vox(:)*dij.ctGrid.resolution.y-currBeam.isoCenter(2);
            zCoordsV       = this.zCoordsV_vox(:)*dij.ctGrid.resolution.z-currBeam.isoCenter(3);
            coordsV        = [xCoordsV yCoordsV zCoordsV];

            xCoordsVdoseGrid = this.xCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.x-currBeam.isoCenter(1);
            yCoordsVdoseGrid = this.yCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.y-currBeam.isoCenter(2);
            zCoordsVdoseGrid = this.zCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.z-currBeam.isoCenter(3);
            coordsVdoseGrid  = [xCoordsVdoseGrid yCoordsVdoseGrid zCoordsVdoseGrid];

            % Get Rotation Matrix
            % Do not transpose matrix since we usage of row vectors &
            % transformation of the coordinate system need double transpose

            currBeam.rotMat_system_T = matRad_getRotationMatrix(currBeam.gantryAngle,currBeam.couchAngle);

            % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
            rot_coordsV         = coordsV*currBeam.rotMat_system_T;
            rot_coordsVdoseGrid = coordsVdoseGrid*currBeam.rotMat_system_T;
            
            rot_coordsV = rot_coordsV - currBeam.sourcePoint_bev;
            rot_coordsVdoseGrid = rot_coordsVdoseGrid - currBeam.sourcePoint_bev;

            % calculate geometric distances
            geoDistVdoseGrid(1:ct.numOfCtScen)= {sqrt(sum(rot_coordsVdoseGrid.^2,2))};

            % Calculate radiological depth cube
            matRad_cfg.dispInfo('matRad: calculate radiological depth cube... ');

            ct.cube = this.cubeWED;
            if this.keepRadDepthCubes
                [radDepthVctGrid, currBeam.radDepthCube] = matRad_rayTracing(currBeam,ct,this.VctGrid,rot_coordsV,this.effectiveLateralCutOff);

                currBeam.radDepthCube = cellfun(@(rD) matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y,   dij.ctGrid.z, rD, ...
                    dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest'),currBeam.radDepthCube,'UniformOutput',false);
                this.radDepthCubes(i,:) = currBeam.radDepthCube(:);
            else
                radDepthVctGrid = matRad_rayTracing(currBeam,ct,this.VctGrid,rot_coordsV,this.effectiveLateralCutOff);
            end
            
            % interpolate radiological depth cube to dose grid resolution
            radDepthVdoseGrid = this.interpRadDepth(ct,1:ct.numOfCtScen,this.VctGrid,this.VdoseGrid,dij.ctGrid,dij.doseGrid,radDepthVctGrid);
            
            % limit rotated coordinates to positions where ray tracing is availabe
            %radDepthsMat = cellfun(@(radDepthCube) matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y,   dij.ctGrid.z,radDepthCube,dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest'),radDepthsMat,'UniformOutput',false);
            
            %Find valid coordinates
            coordIsValid = cellfun(@isfinite, radDepthVdoseGrid,'UniformOutput',false); %Reduce coordinates for finite values
            currBeam.validCoords = cellfun(@and,coordIsValid,this.VdoseGridScenIx,'UniformOutput',false); %Reduce coordinates according to scenario
            currBeam.validCoordsAll = any(cell2mat(coordIsValid),2);
            
            currBeam.radDepths = radDepthVdoseGrid;
            currBeam.geoDepths = geoDistVdoseGrid;
            currBeam.bevCoords = rot_coordsVdoseGrid;      

            % compute SSDs
            currBeam = matRad_computeSSD(currBeam,ct,'densityThreshold',this.ssdDensityThreshold);

            matRad_cfg.dispInfo('done.\n');
           
            %Reinitialize Progress:
            %matRad_progress(1,1000);
        end
        
        function radDepthVdoseGrid = interpRadDepth(~,ct,ctScen,V,Vcoarse,ctGrid,doseGrid,radDepthVctGrid)                        
            for i = 1:numel(ctScen)
                ctScenNum = ctScen(i);

                radDepthCube                = NaN*ones(ct.cubeDim);
                radDepthCube(V(~isnan(radDepthVctGrid{1}))) = radDepthVctGrid{ctScenNum}(~isnan(radDepthVctGrid{1}));

                % interpolate cube - cube is now stored in Y X Z
                coarseRadDepthCube          = matRad_interp3(ctGrid.x,ctGrid.y',ctGrid.z,radDepthCube,doseGrid.x,doseGrid.y',doseGrid.z);
                radDepthVdoseGrid{ctScenNum}  = coarseRadDepthCube(Vcoarse);
            end
        end
        
        function ray = initRay(this,currBeam,j)
            ray = currBeam.ray(j);

            ray.beamIndex = currBeam.beamIndex;
            ray.rayIndex  = j;
            ray.isoCenter = currBeam.isoCenter;
           
            if ~isfield(currBeam,'numOfBixelsPerRay')
                ray.numOfBixels = 1;
            else
                ray.numOfBixels = currBeam.numOfBixelsPerRay(j);
            end

            ray.sourcePoint_bev = currBeam.sourcePoint_bev;
            ray.SAD = currBeam.SAD;
            ray.bixelWidth = currBeam.bixelWidth;

            ray = this.getRayGeometryFromBeam(ray,currBeam);

            %ray = this.computeRaySSD(this,ray); %Is already done in the
            %initBeam function
        end

        function scenRay = extractSingleScenarioRay(this,ray,scenIdx)             
            
             %Gets number of scenario
             scenNum = this.multScen.scenNum(scenIdx);
             ctScen = this.multScen.linearMask(scenNum,1);
            
             %First, create a ray of the
             %specific scenario to adapt rad
             %depths
             scenRay = ray;
             scenRay.radDepths = scenRay.radDepths{ctScen};
             scenRay.radDepths = (1+this.multScen.relRangeShift(scenNum))*scenRay.radDepths + this.multScen.absRangeShift(scenNum);
             scenRay.radialDist_sq = scenRay.radialDist_sq{ctScen};
             scenRay.ix = scenRay.ix{ctScen};
            
             if this.multScen.absRangeShift(scenNum) < 0
                 %TODO: better way to handle this?
                 scenRay.radDepths(scenRay.radDepths < 0) = 0;
             end
            
             if isfield(scenRay,'geoDepths')
                 scenRay.geoDepths = scenRay.geoDepths{ctScen};
             end
            
             if isfield(scenRay,'latDists')
                 scenRay.latDists = scenRay.latDists{ctScen};
             end
            
             if isfield(scenRay,'isoLatDists')
                 scenRay.isoLatDists = scenRay.isoLatDists{ctScen};
             end            
        end
        
        function ray = getRayGeometryFromBeam(this,ray,currBeam)
            lateralRayCutOff = this.getLateralDistanceFromDoseCutOffOnRay(ray);

            % Ray tracing for beam i and ray j
            [ix,radialDist_sq,latDists,isoLatDists] = this.calcGeoDists(currBeam.bevCoords, ...
                ray.sourcePoint_bev, ...
                ray.targetPoint_bev, ...
                ray.SAD, ...
                currBeam.validCoordsAll, ...
                lateralRayCutOff);
            
            %Subindex given the relevant indices from the geometric
            %distance calculation
            ray.validCoords = cellfun(@(beamIx) beamIx & ix,currBeam.validCoords,'UniformOutput',false);
            ray.ix = cellfun(@(ixInGrid) this.VdoseGrid(ixInGrid),ray.validCoords,'UniformOutput',false);
            
            %subCoords = cellfun(@(beamIx) beamIx(ix),currBeam.validCoords,'UniformOutput',false);
            %ray.radialDist_sq = cellfun(@(subix) radialDist_sq(subix),radialDist_sq,subCoords);
            ray.radialDist_sq = cellfun(@(beamIx) radialDist_sq(beamIx(ix)),currBeam.validCoords,'UniformOutput',false);
            ray.latDists = cellfun(@(beamIx) latDists(beamIx(ix),:),currBeam.validCoords,'UniformOutput',false);
            ray.isoLatDists = cellfun(@(beamIx) isoLatDists(beamIx(ix),:),currBeam.validCoords,'UniformOutput',false);

            ray.validCoordsAll = any(cell2mat(ray.validCoords),2);
            
            ray.geoDepths = cellfun(@(rD,ix) rD(ix),currBeam.geoDepths,ray.validCoords,'UniformOutput',false); %usually not needed for particle beams
            ray.radDepths = cellfun(@(rD,ix) rD(ix),currBeam.radDepths,ray.validCoords,'UniformOutput',false);
            %ray.ix = currBeam.ixRadDepths(ix);
            %ray.subIxVdoseGrid = currBeam.subIxVdoseGrid(ix);
        end
        
        function lateralRayCutOff = getLateralDistanceFromDoseCutOffOnRay(this,ray)
            lateralRayCutOff = this.effectiveLateralCutOff;
        end         
        
        function dij = fillDij(this,bixel,dij,stf,scenIdx,currBeamIdx,currRayIdx,currBixelIdx,counter)
            % method for filling the dij struct with the computed dose cube
            % last step in bixel dose calculation

            %Only fill if we actually had bixel (indices) to compute
            if ~isempty(bixel) || ~isempty(bixel.ix)
                % Store in temporary containers to limit matrix filling
                names = fieldnames(this.tmpMatrixContainers);
                bixelContainerColIx = mod(counter-1,this.numOfBixelsContainer)+1;
                subScenIdx = cell(ndims(this.multScen.scenMask),1);
                [subScenIdx{:}] = ind2sub(size(this.multScen.scenMask),scenIdx);
                for q = 1:numel(names)
                    qName = names{q};
                    if this.calcDoseDirect                        
                        %We can omit the resetting to zero because we will
                        %use only the indices we write into
                        %this.tmpMatrixContainers.(qName){bixelContainerColIx,1} = zeros(dij.doseGrid.numOfVoxels,1);                        
                        %this.tmpMatrixContainers.(qName){bixelContainerColIx,1}(this.VdoseGrid(bixel.ix)) = bixel.(qName);
                    else
                        this.tmpMatrixContainers.(qName){bixelContainerColIx,subScenIdx{:}} = sparse(bixel.ix,1,bixel.(qName),dij.doseGrid.numOfVoxels,1);
                    end
                end
                
                % Check if we write to the matrix
                if mod(counter,this.numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
                    if ~this.calcDoseDirect
                        if nargin < 8
                            matRad_cfg = MatRad_Config.instance();
                            matRad_cfg.dispError('Total bixel counter not provided in fillDij');
                        end
        
                        dijColIx = (ceil(counter/this.numOfBixelsContainer)-1)*this.numOfBixelsContainer+1:counter;
                        containerIx = 1:bixelContainerColIx;
                        weight = 1;
                    else
                        dijColIx = currBeamIdx;
                        containerIx = 1;
                        if isfield(stf(currBeamIdx).ray(currRayIdx),'weight') && numel(stf(currBeamIdx).ray(currRayIdx).weight) >= currBixelIdx
                            weight = stf(currBeamIdx).ray(currRayIdx).weight(currBixelIdx);
                        else
                            matRad_cfg = MatRad_Config.instance();
                            matRad_cfg.dispError('No weight available for beam %d, ray %d, bixel %d',currBeamIdx,currRayIdx,currBixelIdx);
                        end
                    end
                    
                    % Iterate through all quantities
                    for q = 1:numel(names)
                        qName = names{q};
                        if ~this.calcDoseDirect
                            dij.(qName){scenIdx}(:,dijColIx) = [this.tmpMatrixContainers.(qName){containerIx,subScenIdx{:}}];
                            %Clean container
                            this.tmpMatrixContainers.(qName)(containerIx,subScenIdx{:}) = cell(numel(containerIx,subScenIdx{:}));
                        else
                            %dij.(qName){1}(this.VdoseGrid(bixel.ix),dijColIx) = dij.(qName){1}(this.VdoseGrid(bixel.ix),dijColIx) + weight * this.tmpMatrixContainers.(qName){containerIx,1}(this.VdoseGrid(bixel.ix));
                            dij.(qName){scenIdx}(bixel.ix,dijColIx) = dij.(qName){scenIdx}(bixel.ix,dijColIx) + weight * bixel.(qName);
                        end
                    end
                end
            end

            %Bookkeeping of bixel numbers
            % remember beam and bixel number
            if this.calcDoseDirect
                dij.beamNum(currBeamIdx)    = currBeamIdx;
                dij.rayNum(currBeamIdx)     = currBeamIdx;
                dij.bixelNum(currBeamIdx)   = currBeamIdx;
            else
                dij.beamNum(counter)    = currBeamIdx;
                dij.rayNum(counter)     = currRayIdx;
                dij.bixelNum(counter)   = currBixelIdx;
            end
        end
        
        %{
        function ray = computeRaySSD(this,ray)
            [alpha,~,rho,d12,~] = matRad_siddonRayTracer(ray.isoCenter, ...
                                 ct.resolution, ...
                                 ray.sourcePoint, ...
                                 ray.targetPoint, ...
                                 this.cubeWED(1));
            ixSSD = find(rho{1} > this.ssdDensityThreshold,1,'first');

            
            if isempty(ixSSD)
                matRad_cfg.dispError('ray does not hit patient. Trying to fix afterwards...');
                boolShowWarning = false;
            elseif ixSSD(1) == 1
                matRad_cfg.dispWarning('Surface for SSD calculation starts directly in first voxel of CT!');
                boolShowWarning = false;
            end
            
            % calculate SSD
            ray.SSD = double(d12* alpha(ixSSD));             
        end
        %}

        function dij = finalizeDose(this,dij)
            %TODO: We could also do this by default for all engines, but
            %this would require to add some additional quantity management
            %(i.e., which quantity an engine could compute)
            % remove dose influence for voxels outside of segmentations for every ct
            % scenario
            for i = 1:this.multScen.numOfCtScen
                % generate index set to erase
                ix = setdiff(1:dij.doseGrid.numOfVoxels,this.VdoseGrid);

                for j = 1:this.multScen.totNumShiftScen
                    for k = 1:this.multScen.totNumRangeScen

                        if this.multScen.scenMask(i,j,k)
                            
                            %loop over all used quantities
                            qNames = fieldnames(this.tmpMatrixContainers);
                            for qIx = 1:numel(qNames)
                                dij.(qNames{qIx}){i,j,k}(ix,:) = 0;
                            end
                        end

                    end
                end
            end

            dij = this.finalizeDose@DoseEngines.matRad_DoseEngineBase(dij);
        end
    end

    methods (Static)

        function [ix,rad_distancesSq,latDists,isoLatDists] = ...
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
                rot_coords_temp = rot_coords_bev(radDepthIx,:);
            else
                % Define rotation matrix
                ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
                R   = eye(3) + ssc(cross(a,b)) + ssc(cross(a,b))^2*(1-dot(a,b))/(norm(cross(a,b))^2);

                % Rotate every CT voxel
                rot_coords_temp = rot_coords_bev(radDepthIx,:)*R;
            end

            % Put [0 0 0] position CT in center of the beamlet.
            %latDistsX = rot_coords_temp(:,1) + sourcePoint_bev(1);
            %latDistsZ = rot_coords_temp(:,3) + sourcePoint_bev(3);
            latDists = rot_coords_temp(:,[1 3]) + sourcePoint_bev([1 3]);

            % check of radial distance exceeds lateral cutoff (projected to iso center)
            %rad_distancesSq = latDistsX.^2 + latDistsZ.^2;
            %subsetMask = rad_distancesSq <= (lateralCutOff/SAD)^2 * rot_coords_temp(:,2).^2;           
            
            rad_distancesSq = sum(latDists.^2,2);
            subsetMask = rad_distancesSq <= (lateralCutOff/SAD)^2 * rot_coords_temp(:,2).^2;
            
            %Apply mask for return quantities

            % index list within considered voxels
            %ix = radDepthIx(subsetMask);
            ix = radDepthIx;
            ix(ix) = subsetMask;

            % return radial distances squared
            if nargout > 1
                rad_distancesSq = rad_distancesSq(subsetMask);
            end
            
            %lateral distances in X & Z
            if nargout > 2
                latDists = latDists(subsetMask,:);
            end
            
            % lateral distances projected onto isocenter
            if nargout > 3
                isoLatDists = latDists./rot_coords_temp(subsetMask,2)*SAD;           
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

