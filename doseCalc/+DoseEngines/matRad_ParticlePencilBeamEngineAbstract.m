classdef (Abstract) matRad_ParticlePencilBeamEngineAbstract < DoseEngines.matRad_PencilBeamEngineAbstract
    % matRad_DoseEngineParticlePB:
    %   Implements an engine for particle based dose calculation
    %   For detailed information see superclass matRad_DoseEngine
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2022 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % help edit

    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (SetAccess = public, GetAccess = public)

        calcLET = true;                 % Boolean which defines if LET should be calculated
        calcBioDose = false;            % Boolean which defines if biological dose calculation shoudl be performed (alpha*dose and sqrt(beta)*dose)

        visBoolLateralCutOff = false;   % Boolean switch for visualization during+ LeteralCutOff calculation

        lateralModel = 'auto';          % Lateral Model used. 'auto' uses the most accurate model available (i.e. multiple Gaussians). 'single','double','multi' try to force a singleGaussian or doubleGaussian model, if available
    end

    properties (SetAccess = protected, GetAccess = public)
        constantRBE = NaN;              % constant RBE value

        vTissueIndex;                   %Stores tissue indices available in the matRad base data
    end

    methods
        function this = matRad_ParticlePencilBeamEngineAbstract(pln)
            this = this@DoseEngines.matRad_PencilBeamEngineAbstract(pln);

            % check if bio optimization is needed and set the
            % coresponding boolean accordingly
            % TODO:
            % This should not be handled here as an optimization property
            % We should rather make optimization dependent on what we have
            % decided to calculate here.
            if nargin > 0 
                if (isfield(pln,'propOpt')&& isfield(pln.propOpt,'bioOptimization')&& ...
                    (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') ||...
                    isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) && ...
                    strcmp(pln.radiationMode,'carbon'))
                this.calcBioDose = true;
                elseif strcmp(pln.radiationMode,'protons') && isfield(pln,'propOpt') && isfield(pln.propOpt,'bioOptimization') && isequal(pln.propOpt.bioOptimization,'const_RBExD')
                    this.constantRBE = 1.1;                    
                end
            end
        end
        

    end

    % Should be abstract methods but in order to satisfy the compatibility
    % with OCTAVE we can't use abstract methods. If OCTAVE at some point
    % in the far future implements this feature this should be abstract again.
    methods (Access = protected) %Abstract
        function bixel = calcParticleBixel(this,bixel)
            error('Function needs to be implemented!');
        end
    end

    methods (Access = protected)
        function chooseLateralModel(this)
            fValidateMulti = @(bd) isfield(bd,'sigmaMulti') && isfield(bd,'weightMulti') && ~isempty(bd.sigmaMulti) && ~isempty(bd.weightMulti);
            fValidateDouble = @(bd) isfield(bd,'sigma1') && isfield(bd,'sigma2') && isfield(bd,'weight') && ~isempty(bd.sigma1) && ~isempty(bd.sigma2) && ~isempty(bd.weight);
            fValidateSingle = @(bd) isfield(bd,'sigma') && ~isempty(bd.sigma);
            
            matRad_cfg = MatRad_Config.instance();

            switch this.lateralModel
                case 'single'
                    if ~all(arrayfun(fValidateSingle,this.machine.data))
                        matRad_cfg.dispWarning('Chosen Machine does not support a singleGaussian Pencil-Beam model!');
                        this.lateralModel = 'auto';
                    end
                case 'double'
                    if ~all(arrayfun(fValidateDouble,this.machine.data))
                        matRad_cfg.dispWarning('Chosen Machine does not support a doubleGaussian Pencil-Beam model!');
                        this.lateralModel = 'auto';
                    end
                case 'multi'
                    if ~all(arrayfun(fValidateMulti,this.machine.data))
                        matRad_cfg.dispWarning('Chosen Machine does not support a multiGaussian Pencil-Beam model!');
                        this.lateralModel = 'auto';
                    end
                case 'auto'
                    %Do nothing, will be handled below
                otherwise
                    matRad_cfg.dispError('Lateral model ''%s'' not known!',this.lateralModel);
            end

            %Now check if we need tho chose the lateral model because it
            %was set to auto
            if strcmp(this.lateralModel,'auto') 
                if all(arrayfun(fValidateMulti,this.machine.data))
                    this.lateralModel = 'multi';
                elseif all(arrayfun(fValidateDouble,this.machine.data))
                    this.lateralModel = 'double';
                elseif all(arrayfun(fValidateSingle,this.machine.data))
                    this.lateralModel = 'single';
                else
                    matRad_cfg.dispError('Invalid kernel model!');
                end
            end

            matRad_cfg.dispInfo('Using a %s Gaussian pencil-beam kernel model!\n');
        end
        
        function bixel = computeBixel(this,currRay,k)
            %Initialize Bixel Geometry
            bixel = this.initBixel(currRay,k);
            
            %Compute Bixel
            bixel = this.calcParticleBixel(bixel);
        end

        function bixel = initBixel(this,currRay,k)
            % matRad initialize general bixel geometry for particle dose calc
            %
            % call
            %   bixel = this.initBixel(currRay,k)
            
            bixel = struct();
            bixel.beamIndex = currRay.beamIndex;
            bixel.rayIndex = currRay.rayIndex;
            bixel.bixelIndex = k;
                     
            %First we get metadata: MU, corresponding base data entry, etc.
            % extract MU data if present (checks for downwards compatability)
            bixel.minMU = 0;
            if isfield(currRay,'minMU')
                bixel.minMU = currRay.minMU(k);
            end

            bixel.maxMU = Inf;
            if isfield(currRay,'maxMU')
                bixel.maxMU = currRay.maxMU(k);
            end

            bixel.numParticlesPerMU = 1e6;
            if isfield(currRay,'numParticlesPerMU')
                bixel.numParticlesPerMU = currRay.numParticlesPerMU(k);
            end

            % find energy index in base data
            energyIx = find(this.round2(currRay.energy(k),4) == this.round2([this.machine.data.energy],4));
            bixel.energyIx = energyIx;
            bixel.baseData = this.machine.data(energyIx);
            
            bixel.rangeShifter = currRay.rangeShifter(k);
            
            bixel.SSD = currRay.SSD;

            % Compute initial spotWidth
            bixel.sigmaIniSq = currRay.sigmaIni(k).^2;

            % Apply beam modifiers
            bixel = this.getBeamModifiers(bixel);
           
            %Gets bixel.ix (voxel indices) and bixel.subIx (logical
            %indices to be used)  after cutoff. Storing these allows us to
            %use indexing for performance and avoid to many copies
            bixel = this.getBixelIndicesOnRay(bixel,currRay);

            if isempty(bixel.ix)
                return;
            end     

            % Get quantities 1:1 from ray. Here we trust Matlab's memory
            % management to not copy the arrays until they are modified.
            % This allows us to efficiently access them by indexing in the
            % bixel computation
            bixel.radialDist_sq = currRay.radialDist_sq(bixel.subRayIx);
            bixel.radDepths = currRay.radDepths(bixel.subRayIx);
            if this.calcBioDose
                bixel.vTissueIndex = currRay.vTissueIndex(bixel.subRayIx);
            end
        
        end

        function X = interpolateKernelsInDepth(this,bixel)
            baseData = bixel.baseData;
            
            depths = baseData.depths;

            % add potential offset
            if isfield(baseData,'offset')
                depths = depths + baseData.offset - bixel.radDepthOffset;
            end
            
            % calculate particle dose for bixel k on ray j of beam i
            % convert from MeV cm^2/g per primary to Gy mm^2 per 1e6 primaries
            conversionFactor = 1.6021766208e-02;
            
            %Find all values we need to interpolate
            X.Z = conversionFactor*baseData.Z;
            
            %Lateral Kernel Model
            switch this.lateralModel
                case 'single'
                    X.sigma  = baseData.sigma;
                case 'double'
                    X.sigma1 = baseData.sigma1;
                    X.sigma2 = baseData.sigma2;
                    X.weight = baseData.weight;
                case 'multi'
                    X.weightMulti = baseData.weightMulti;
                    X.sigmaMulti = baseData.sigmaMulti;
                otherwise
                    %Sanity check
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid Lateral Model');
            end
                   
            % LET
            if this.calcLET
                X.LET = baseData.LET;
            end

            % bioDose
            if this.calcBioDose
                X.alpha = baseData.alpha;
                X.beta = baseData.beta;
            end   
            
            X = structfun(@(v) matRad_interp1(depths,v,bixel.radDepths,'nearest'),X,'UniformOutput',false); %Extrapolate to zero?           
        end

        % We override this function to boost efficiency a bit (latDistX & Z
        % not needed)
        function ray = getRayGeometryFromBeam(this,ray,currBeam)
            lateralRayCutOff = this.getLateralDistanceFromDoseCutOffOnRay(ray);

            % Ray tracing for beam i and ray j
            [ray.ix,ray.radialDist_sq] = this.calcGeoDists(currBeam.rot_coordsVdoseGrid, ...
                ray.sourcePoint_bev, ...
                ray.targetPoint_bev, ...
                ray.SAD, ...
                currBeam.ixRadDepths, ...
                lateralRayCutOff);
            
            ray.radDepths = currBeam.radDepthVdoseGrid{1}(ray.ix);
        end

        function [currBixel] = getBixelIndicesOnRay(this,currBixel,currRay)
            
            % create offset vector to account for additional offsets modelled in the base data and a potential
            % range shifter. In the following, we only perform dose calculation for voxels having a radiological depth
            % that is within the limits of the base data set (-> machine.data(i).dephts). By this means, we only allow
            % interpolations in this.calcParticleDoseBixel() and avoid extrapolations.
            %urrBixel.offsetRadDepth = currBixel.baseData.offset + currBixel.radDepthOffset;
            tmpOffset = currBixel.baseData.offset - currBixel.radDepthOffset;

            % find depth depended lateral cut off
            if this.dosimetricLateralCutOff == 1
                currIx = currRay.radDepths <= currBixel.baseData.depths(end) + tmpOffset;
            elseif this.dosimetricLateralCutOff < 1 && this.dosimetricLateralCutOff > 0
                
                %{
                % Old two step clipping. Seems to be slower than doing it
                at once using nearest neighbor extrapolation
                % perform rough 2D clipping
                currIx = currRay.radDepths <= currBixel.baseData.depths(end) + tmpOffset & ...
                    currRay.radialDist_sq <= max(currBixel.baseData.LatCutOff.CutOff.^2);

                % peform fine 2D clipping
                if length(currBixel.baseData.LatCutOff.CutOff) > 1
                    currIx(currIx) = matRad_interp1((currBixel.baseData.LatCutOff.depths + tmpOffset)',...
                        (currBixel.baseData.LatCutOff.CutOff.^2)', currRay.radDepths(currIx)) >= currRay.radialDist_sq(currIx);
                end
                %}                

                if length(currBixel.baseData.LatCutOff.CutOff) > 1
                    %currIx = matRad_interp1((currBixel.baseData.LatCutOff.depths + tmpOffset)',(currBixel.baseData.LatCutOff.CutOff.^2)', currRay.radDepths,'nearest') >= currRay.radialDist_sq & currRay.radDepths <= currBixel.baseData.depths(end) + tmpOffset;
                    currIx = matRad_interp1((currBixel.baseData.LatCutOff.depths + tmpOffset)',(currBixel.baseData.LatCutOff.CutOff.^2)', currRay.radDepths) >= currRay.radialDist_sq & currRay.radDepths <= currBixel.baseData.depths(end) + tmpOffset;
                else
                    currIx = currBixel.baseData.LatCutOff.CutOff.^2 >= currRay.radialDist_sq & currRay.radDepths <= currBixel.baseData.depths(end) + tmpOffset;
                end
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Cutoff must be a value between 0 and 1!')
            end

            currBixel.subRayIx = currIx;
            currBixel.ix = currRay.ix(currIx);
        end

        function currBixel = getBeamModifiers(this,currBixel)
            radDepthOffset = 0;
            addSigmaSq = 0;
            % consider range shifter for protons if applicable
            if currBixel.rangeShifter.eqThickness > 0
                %TODO: We should do this check in dose calc initialization
                %instead to spare some time
                if ~strcmp(this.machine.meta.radiationMode,'protons')
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispWarning('Range shifter not valid for irradiation with particle other than protons!');
                end

                % compute!
                sigmaRashi = matRad_calcSigmaRashi(currBixel.baseData.energy, ...
                    currBixel.rangeShifter, ...
                    currBixel.SSD);

                % add to initial sigma in quadrature
                addSigmaSq = addSigmaSq + sigmaRashi.^2;
                radDepthOffset = radDepthOffset + currBixel.rangeShifter.eqThickness;              

            end

            currBixel.addSigmaSq = addSigmaSq;
            currBixel.radDepthOffset = radDepthOffset;
        end
        
        function [dij,ct,cst,stf] = calcDoseInit(this,ct,cst,stf)
            % modified inherited method of the superclass DoseEngine,
            % containing intialization which are specificly needed for
            % pencil beam calculation and not for other engines

            [dij,ct,cst,stf] = calcDoseInit@DoseEngines.matRad_PencilBeamEngineAbstract(this,ct,cst,stf);
            
            %Choose the lateral pencil beam model
            this.chooseLateralModel();

            if ~isnan(this.constantRBE)
                dij.RBE = this.constantRBE;
            end

            % Load biologicla base data if needed
            if this.calcBioDose
                dij = this.loadBiologicalBaseData(cst,dij);
                % allocate alpha and beta dose container and sparse matrices in the dij struct,
                % for more informations see corresponding method
                dij = this.allocateBioDoseContainer(dij);
            end           

            % allocate LET containner and let sparse matrix in dij struct
            if this.calcLET
                if isfield(this.machine.data,'LET')
                    dij = this.allocateLETContainer(dij);
                else
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispWarning('LET not available and will not be computed!');
                    this.calcLET = false;
                end
            end

            % lateral cutoff for raytracing and geo calculations
            this.effectiveLateralCutOff = this.geometricLateralCutOff;
        end

        function currBeam = initBeam(this,dij,ct,cst,stf,i)
            currBeam = initBeam@DoseEngines.matRad_PencilBeamEngineAbstract(this,dij,ct,cst,stf,i);
            
            %currBeam.rot_coordsVdoseGrid = currBeam.rot_coordsVdoseGrid(~isnan(currBeam.radDepthVdoseGrid{1}),:);
            maxEnergy = max([currBeam.ray.energy]);
            maxEnergyIx = find(this.round2(maxEnergy,4) == this.round2([this.machine.data.energy],4));
            raShis = [currBeam.ray.rangeShifter];
            minRaShi = min([raShis.eqThickness]);
            radDepthOffset = this.machine.data(maxEnergyIx).offset + minRaShi;

            % limit rotated coordinates to positions where ray tracing is availabe
            currBeam.ixRadDepths = find(currBeam.radDepthVdoseGrid{1} < (this.machine.data(maxEnergyIx).depths(end) - radDepthOffset));% & ~isnan(currBeam.radDepthVdoseGrid{1}));
            
            %currBeam.radDepthVdoseGrid{1} = currBeam.radDepthVdoseGrid{1}(currBeam.ixRadDepths);
            currBeam.rot_coordsVdoseGrid = currBeam.rot_coordsVdoseGrid(currBeam.ixRadDepths,:);

            %Precompute CutOff
            this.calcLateralParticleCutOff(this.dosimetricLateralCutOff,stf(i));
        end
        
        function dij = loadBiologicalBaseData(this,cst,dij)
            matRad_cfg = MatRad_Config.instance();
            if isfield(this.machine.data,'alphaX') && isfield(this.machine.data,'betaX')
                matRad_cfg.dispInfo('matRad: loading biological base data... ');
                this.vTissueIndex    = zeros(size(this.VdoseGrid,1),1);
                dij.ax              = zeros(dij.doseGrid.numOfVoxels,1);
                dij.bx              = zeros(dij.doseGrid.numOfVoxels,1);

                cstDownsampled = matRad_setOverlapPriorities(cst);

                % resizing cst to dose cube resolution
                cstDownsampled = matRad_resizeCstToGrid(cstDownsampled,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,...
                    dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);
                % retrieve photon LQM parameter for the current dose grid voxels
                [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cstDownsampled,dij.doseGrid.numOfVoxels,1,this.VdoseGrid);

                for i = 1:size(cstDownsampled,1)

                    % check if cst is compatiable
                    if ~isempty(cstDownsampled{i,5}) && isfield(cstDownsampled{i,5},'alphaX') && isfield(cstDownsampled{i,5},'betaX')

                        % check if base data contains alphaX and betaX
                        IdxTissue = find(ismember(this.machine.data(1).alphaX,cstDownsampled{i,5}.alphaX) & ...
                            ismember(this.machine.data(1).betaX,cstDownsampled{i,5}.betaX));

                        % check consitency of biological baseData and cst settings
                        if ~isempty(IdxTissue)
                            isInVdoseGrid = ismember(this.VdoseGrid,cstDownsampled{i,4}{1});
                            this.vTissueIndex(isInVdoseGrid) = IdxTissue;
                        else
                            matRad_cfg.dispError('Biological base data and cst are inconsistent!');
                        end

                    else
                        this.vTissueIndex(row) = 1;
                        matRad_cfg.dispInfo('Tissue type of %s was set to 1\n',cstDownsampled{i,2});
                    end
                end
                dij.vTissueIndex = this.vTissueIndex;
                matRad_cfg.dispInfo('done.\n');
            else
                matRad_cfg.dispError('Base data is missing alphaX and/or betaX!');
            end
        end

        function dij = allocateBioDoseContainer(this,dij)
            % allocate space for container used in bio optimization
            dij = this.allocateQuantityMatrixContainers(dij,{'mAlphaDose','mSqrtBetaDose'});
        end

        function dij = allocateLETContainer(this,dij)
            % allocate space for container used in LET calculation

            % get MatLab Config instance for displaying warings
            matRad_cfg = MatRad_Config.instance();
            if isfield(this.machine.data,'LET')
                dij = this.allocateQuantityMatrixContainers(dij,{'mLETDose'});         
            else
                matRad_cfg.dispWarning('LET not available in the machine data. LET will not be calculated.');
            end

        end
        
        function calcLateralParticleCutOff(this,cutOffLevel,stfElement)
            % matRad function to calculate a depth dependend lateral cutoff
            % for each pristine particle beam
            %
            % call
            %   this.calcLateralParticleCutOff(cutOffLevel,stf)
            %
            % input
            %   this:        current engine object includes machine base data file
            %   cutOffLevel:    cut off level - number between 0 and 1
            %   stfElement:  matRad steering information struct for a single beam
            %
            % output
            %   machine:    	changes in the object property machine base data file including an additional field representing the lateral
            %                    cutoff
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
            
            matRad_cfg = MatRad_Config.instance();
            
            %Sanity Checks
            if numel(stfElement) > 1
                matRad_cfg.dispError('CutOff can only be precalculated for a single element, but you provided steering information for multiple beams!');
            end

            if cutOffLevel <= 0.98
                matRad_cfg.dispWarning('a lateral cut off below 0.98 may result in an inaccurate dose calculation');
            end

            if (cutOffLevel < 0 || cutOffLevel > 1)
                matRad_cfg.dispWarning('lateral cutoff is out of range - using default cut off of 0.99')
                cutOffLevel = 0.99;
            end


            matRad_cfg.dispInfo('matRad: calculate lateral cutoff...');           
            conversionFactor = 1.6021766208e-02;

            % function handle for calculating depth dose for APM
            sumGauss = @(x,mu,SqSigma,w) ((1./sqrt(2*pi*ones(numel(x),1) * SqSigma') .* ...
                exp(-bsxfun(@minus,x,mu').^2 ./ (2* ones(numel(x),1) * SqSigma' ))) * w);

            % define some variables needed for the cutoff calculation
            vX = [0 logspace(-1,3,1200)]; % [mm]

            % integration steps
            r_mid          = 0.5*(vX(1:end-1) +  vX(2:end))'; % [mm]
            dr             = (vX(2:end) - vX(1:end-1))';
            radialDist_sq  = r_mid.^2;

            % number of depth points for which a lateral cutoff is determined
            numDepthVal    = 35;

            % helper function for energy selection
            round2 = @(a,b)round(a*10^b)/10^b;

            % extract SSD for each bixel
            vSSD = ones(1,length([stfElement.ray(:).energy]));
            cnt = 1;
            for i  = 1:length(stfElement.ray)
                vSSD(cnt:cnt+numel([stfElement.ray(i).energy])-1) = stfElement.ray(i).SSD;
                cnt = cnt + numel(stfElement.ray(i).energy);
            end

            % setup energy, focus index, sigma look up table - only consider unique rows
            [energySigmaLUT,ixUnique]  = unique([[stfElement.ray(:).energy]; [stfElement.ray(:).focusIx] ; vSSD]','rows');
            rangeShifterLUT = [stfElement.ray(:).rangeShifter];
            rangeShifterLUT = rangeShifterLUT(1,ixUnique);

            % find the largest inital beam width considering focus index, SSD and range shifter for each individual energy
            for i = 1:size(energySigmaLUT,1)

                % find index of maximum used energy (round to keV for numerical reasons
                energyIx = max(round2(energySigmaLUT(i,1),4)) == round2([this.machine.data.energy],4);

                currFoci = energySigmaLUT(i,2);
                sigmaIni = matRad_interp1(this.machine.data(energyIx).initFocus.dist(currFoci,:)',...
                    this.machine.data(energyIx).initFocus.sigma(currFoci,:)',...
                    energySigmaLUT(i,3));
                sigmaIni_sq = sigmaIni^2;

                % consider range shifter for protons if applicable
                if  strcmp(this.machine.meta.radiationMode,'protons') && rangeShifterLUT(i).eqThickness > 0  && ~strcmp(this.machine.meta.machine,'Generic')

                    %get max range shift
                    sigmaRashi = matRad_calcSigmaRashi(this.machine.data(energyIx).energy, ...
                        rangeShifterLUT(i), ...
                        energySigmaLUT(i,3));

                    % add to initial sigma in quadrature
                    sigmaIni_sq = sigmaIni_sq +  sigmaRashi.^2;

                end

                energySigmaLUT(i,4) = sigmaIni_sq;

            end

            % find for each individual energy the broadest inital beam width
            uniqueEnergies                = unique(energySigmaLUT(:,1));
            largestSigmaSq4uniqueEnergies = NaN * ones(numel(uniqueEnergies),1);
            ix_Max                        = NaN * ones(numel(uniqueEnergies),1);
            for i = 1:numel(uniqueEnergies)
                [largestSigmaSq4uniqueEnergies(i), ix_Max(i)] = max(energySigmaLUT(uniqueEnergies(i) == energySigmaLUT(:,1),4));
            end

            % get energy indices for looping
            vEnergiesIx = find(ismember([this.machine.data(:).energy],uniqueEnergies(:,1)));
            cnt         = 0;

            % loop over all entries in the machine.data struct
            for energyIx = vEnergiesIx

                % set default depth cut off - finite value will be set during first iteration
                depthDoseCutOff = inf;

                % get the current integrated depth dose profile
                if isstruct(this.machine.data(energyIx).Z)
                    idd_org = sumGauss(this.machine.data(energyIx).depths,this.machine.data(energyIx).Z.mean,...
                        this.machine.data(energyIx).Z.width.^2,...
                        this.machine.data(energyIx).Z.weight) * conversionFactor;
                else
                    idd_org = this.machine.data(energyIx).Z * conversionFactor;
                end

                [~,peakIxOrg] = max(idd_org);

                % get indices for which a lateral cutoff should be calculated
                cumIntEnergy = cumtrapz(this.machine.data(energyIx).depths,idd_org);

                peakTailRelation   = 0.5;
                numDepthValToPeak  = ceil(numDepthVal*peakTailRelation);                                                                          % number of depth values from 0 to peak position
                numDepthValTail    = ceil(numDepthVal*(1-peakTailRelation));                                                                      % number of depth values behind peak position
                energyStepsToPeak  = cumIntEnergy(peakIxOrg)/numDepthValToPeak;
                energyStepsTail    = (cumIntEnergy(end)-cumIntEnergy(peakIxOrg))/numDepthValTail;
                % make sure to include 0, peak position and end position
                vEnergySteps       = unique([0:energyStepsToPeak:cumIntEnergy(peakIxOrg) cumIntEnergy(peakIxOrg) ...
                    cumIntEnergy(peakIxOrg+1):energyStepsTail:cumIntEnergy(end) cumIntEnergy(end)]);

                [cumIntEnergy,ix] = unique(cumIntEnergy);
                depthValues       = matRad_interp1(cumIntEnergy,this.machine.data(energyIx).depths(ix),vEnergySteps);

                if isstruct(this.machine.data(energyIx).Z)
                    idd = sumGauss(depthValues,this.machine.data(energyIx).Z.mean,...
                        this.machine.data(energyIx).Z.width.^2,...
                        this.machine.data(energyIx).Z.weight) * conversionFactor;
                else
                    idd  = matRad_interp1(this.machine.data(energyIx).depths,this.machine.data(energyIx).Z,depthValues) * conversionFactor;
                end

                cnt = cnt +1 ;
                % % calculate dose in spot
                baseData                   = this.machine.data(energyIx);
                baseData.LatCutOff.CompFac = 1;

                for j = 1:numel(depthValues)

                    % save depth value
                    this.machine.data(energyIx).LatCutOff.depths(j) = depthValues(j);

                    if cutOffLevel == 1
                        this.machine.data(energyIx).LatCutOff.CompFac   = 1;
                        this.machine.data(energyIx).LatCutOff.CutOff(j) = Inf;
                    else
                        bixel.energyIx = energyIx;
                        bixel.baseData = baseData;
                        bixel.radialDist_sq = radialDist_sq;
                        bixel.sigmaIniSq = largestSigmaSq4uniqueEnergies(cnt);
                        bixel.radDepths = (depthValues(j) + baseData.offset) * ones(size(radialDist_sq));
                        bixel.vTissueIndex = ones(size(bixel.radDepths));
                        bixel.subRayIx = true(size(bixel.radDepths));
                        bixel.ix = find(bixel.subRayIx);
                        bixel.radDepthOffset = 0;
                        bixel.addSigmaSq = 0;
                        
                        % calculate dose
                        bixel = this.calcParticleBixel(bixel);
                        dose_r = bixel.physicalDose;

                        cumArea = cumsum(2*pi.*r_mid.*dose_r.*dr);
                        relativeTolerance = 0.5; %in [%]
                        if abs((cumArea(end)./(idd(j)))-1)*100 > relativeTolerance
                            matRad_cfg.dispWarning('LateralParticleCutOff: shell integration is wrong !')
                        end

                        IX = find(cumArea >= idd(j) * cutOffLevel,1, 'first');
                        this.machine.data(energyIx).LatCutOff.CompFac = cutOffLevel^-1;

                        if isempty(IX)
                            depthDoseCutOff = Inf;
                            matRad_cfg.dispWarning('LateralParticleCutOff: Couldnt find lateral cut off !')
                        elseif isnumeric(IX)
                            depthDoseCutOff = r_mid(IX);
                        end

                        this.machine.data(energyIx).LatCutOff.CutOff(j) = depthDoseCutOff;

                    end
                end
            end

            matRad_cfg.dispInfo('done.\n');

            %% visualization
            if this.visBoolLateralCutOff

                % determine which pencil beam should be plotted
                subIx    = ceil(numel(vEnergiesIx)/2);
                energyIx = vEnergiesIx(subIx);

                baseData       = this.machine.data(energyIx);
                focusIx        = energySigmaLUT(ix_Max(subIx),2);
                maxSSD         = energySigmaLUT(ix_Max(subIx),3);
                rangeShifter   = rangeShifterLUT(ix_Max(subIx));
                TmpCompFac     = baseData.LatCutOff.CompFac;
                baseData.LatCutOff.CompFac = 1;

                % plot 3D cutoff at one specific depth on a rather sparse grid
                sStep         = 0.5;
                vLatX         = -100 : sStep : 100; % [mm]
                dimX          = numel(vLatX);
                midPos        = round(length(vLatX)/2);
                [X,Y]         = meshgrid(vLatX,vLatX);

                radDepths     = [0:sStep:this.machine.data(energyIx).depths(end)] + this.machine.data(energyIx).offset;
                radialDist_sq = (X.^2 + Y.^2);
                radialDist_sq = radialDist_sq(:);
                mDose         = zeros(dimX,dimX,numel(radDepths));
                vDoseInt      = zeros(numel(radDepths),1);

                for kk = 1:numel(radDepths)

                    % calculate initial focus sigma
                    sigmaIni = matRad_interp1(this.machine.data(energyIx).initFocus.dist(focusIx,:)', ...
                        this.machine.data(energyIx).initFocus.sigma(focusIx,:)',maxSSD);
                    sigmaIni_sq = sigmaIni^2;

                    % consider range shifter for protons if applicable
                    if rangeShifter.eqThickness > 0 && strcmp(pln.radiationMode,'protons')

                        % compute!
                        sigmaRashi = matRad_calcSigmaRashi(this.machine.data(energyIx).energy,rangeShifter,maxSSD);

                        % add to initial sigma in quadrature
                        sigmaIni_sq = sigmaIni_sq +  sigmaRashi^2;

                    end
                    
                    bixel.energyIx = energyIx;
                    bixel.baseData = baseData;
                    bixel.radialDist_sq = radialDist_sq;
                    bixel.sigmaIniSq = sigmaIni_sq;
                    bixel.radDepths = radDepths(kk)*ones(size(bixel.radialDist_sq));
                    bixel.vTissueIndex = ones(size(bixel.radDepths));
                    bixel.subRayIx = true(size(bixel.radDepths));
                    bixel.ix = find(bixel.subRayIx);
                    bixel.radDepthOffset = 0;
                    bixel.addSigmaSq = 0;

                    bixel = this.calcParticleBixel(bixel);
                    
                    mDose(:,:,kk) = reshape(bixel.physicalDose,[dimX dimX]);

                    [~,IX]           = min(abs((this.machine.data(energyIx).LatCutOff.depths + this.machine.data(energyIx).offset) - radDepths(kk)));
                    TmpCutOff        = this.machine.data(energyIx).LatCutOff.CutOff(IX);
                    vXCut            = vX(vX<=TmpCutOff);

                    % integration steps
                    r_mid_Cut        = (0.5*(vXCut(1:end-1) +  vXCut(2:end)))'; % [mm]
                    dr_Cut           = (vXCut(2:end) - vXCut(1:end-1))';
                    radialDist_sqCut = r_mid_Cut.^2;
                    
                    bixel.radialDist_sq = radialDist_sqCut(:);
                    bixel.radDepths = radDepths(kk)*ones(size(bixel.radialDist_sq));
                    bixel.vTissueIndex = ones(size(bixel.radDepths));
                    bixel.subRayIx = true(size(bixel.radDepths));

                    bixel = this.calcParticleBixel(bixel);
                    dose_r_Cut       = bixel.physicalDose;

                    cumAreaCut = cumsum(2*pi.*r_mid_Cut.*dose_r_Cut.*dr_Cut);

                    if ~isempty(cumAreaCut)
                        vDoseInt(kk) = cumAreaCut(end);
                    end
                end

                % obtain maximum dose
                if isstruct(this.machine.data(energyIx).Z)
                    idd = sumGauss(depthValues,this.machine.data(energyIx).Z.mean,...
                        this.machine.data(energyIx).Z.width.^2,...
                        this.machine.data(energyIx).Z.weight) * conversionFactor;
                else
                    idd  = matRad_interp1(this.machine.data(energyIx).depths,this.machine.data(energyIx).Z,depthValues) * conversionFactor;
                end

                [~,peakixDepth] = max(idd);

                bixel.energyIx = energyIx;
                bixel.baseData = baseData;
                bixel.radialDist_sq = 0;
                bixel.sigmaIniSq = sigmaIni_sq;
                bixel.radDepths = baseData.depths(peakixDepth);
                bixel.vTissueIndex = ones(size(bixel.radDepths));
                bixel.subIx = true;
                bixel.radDepthOffset = 0;
                bixel.addSigmaSq = 0;
                
                bixel = this.calcParticleBixel(bixel);
                dosePeakPos = bixel.physicalDose;

                vLevelsDose = dosePeakPos.*[0.01 0.05 0.1 0.9];
                doseSlice   = squeeze(mDose(midPos,:,:));
                figure,set(gcf,'Color',[1 1 1]);
                subplot(311),h=imagesc(squeeze(mDose(midPos,:,:)));hold on;
                set(h,'AlphaData', .8*double(doseSlice>0));
                contour(doseSlice,vLevelsDose,'LevelListMode','manual','LineWidth',2);hold on

                ax = gca;
                ax.XTickLabelMode = 'manual';
                ax.XTickLabel     = strsplit(num2str(ax.XTick*sStep + this.machine.data(energyIx).offset),' ')';
                ax.YTickLabelMode = 'manual';
                ax.YTickLabel     = strsplit(num2str(ax.YTick*sStep + this.machine.data(energyIx).offset),' ')';

                plot(1+(this.machine.data(energyIx).LatCutOff.depths)*sStep^-1,...
                    this.machine.data(energyIx).LatCutOff.CutOff * sStep^-1 + midPos,'rx');

                legend({'isodose 1%,5%,10% 90%','calculated cutoff'}) ,colorbar,set(gca,'FontSize',12),xlabel('z [mm]'),ylabel('x [mm]');

                entry = this.machine.data(energyIx);
                if isstruct(entry.Z)
                    idd = sumGauss(entry.depths,entry.Z.mean,entry.Z.width.^2,entry.Z.weight);
                else
                    idd = this.machine.data(energyIx).Z;
                end
                subplot(312),plot(this.machine.data(energyIx).depths,idd*conversionFactor,'k','LineWidth',2),grid on,hold on
                plot(radDepths - this.machine.data(energyIx).offset,vDoseInt,'r--','LineWidth',2),hold on,
                plot(radDepths - this.machine.data(energyIx).offset,vDoseInt * TmpCompFac,'bx','LineWidth',1),hold on,
                legend({'original IDD',['cut off IDD at ' num2str(cutOffLevel) '%'],'cut off IDD with compensation'},'Location','northwest'),
                xlabel('z [mm]'),ylabel('[MeV cm^2 /(g * primary)]'),set(gca,'FontSize',12)

                totEnergy        = trapz(this.machine.data(energyIx).depths,idd*conversionFactor) ;
                totEnergyCutOff  = trapz(radDepths,vDoseInt * TmpCompFac) ;
                relDiff          =  ((totEnergy/totEnergyCutOff)-1)*100;
                title(['rel diff of integral dose ' num2str(relDiff) '%']);
                baseData.LatCutOff.CompFac = TmpCompFac;

                subplot(313),
                if isfield(this.machine.data(energyIx),'sigma1')
                    yyaxis left;
                    plot(this.machine.data(energyIx).LatCutOff.depths,this.machine.data(energyIx).LatCutOff.CutOff,'LineWidth',2),hold on
                    plot(this.machine.data(energyIx).depths,(this.machine.data(energyIx).sigma1),':','LineWidth',2),grid on,hold on,ylabel('mm')
                    yyaxis right;
                    plot(this.machine.data(energyIx).depths,(this.machine.data(energyIx).sigma2),'-.','LineWidth',2),grid on,hold on,ylabel('mm')
                    legend({'Cutoff','sigma1','sigma2'});
                else
                    yyaxis left;plot(this.machine.data(energyIx).LatCutOff.depths,this.machine.data(energyIx).LatCutOff.CutOff,'LineWidth',2),hold on,ylabel('mm')
                    yyaxis right;subplot(313),plot(this.machine.data(energyIx).depths,this.machine.data(energyIx).sigma,'LineWidth',2),grid on,hold on
                    legend({'Cutoff','sigma'});ylabel('mm')
                end

                set(gca,'FontSize',12),xlabel('z [mm]'),  ylabel('mm')

                % plot cutoff of different energies
                figure,set(gcf,'Color',[1 1 1]);
                cnt = 1;
                for i = vEnergiesIx
                    plot(this.machine.data(i).LatCutOff.depths,this.machine.data(i).LatCutOff.CutOff,'LineWidth',1.5),hold on
                    cellLegend{cnt} = [num2str(this.machine.data(i).energy) ' MeV'];
                    cnt = cnt + 1;
                end
                grid on, grid minor,xlabel('depth in [mm]'),ylabel('lateral cutoff in [mm]')
                title(['cutoff level = ' num2str(cutOffLevel)]),
                ylim = get(gca,'Ylim');    set(gca,'Ylim',[0 ylim(2)+3]),    legend(cellLegend)
            end
        end
    
        function ray = initRay(this,beam,j)           
            ray = initRay@DoseEngines.matRad_PencilBeamEngineAbstract(this,beam,j);

            % calculate initial sigma for all bixel on current ray
            ray.sigmaIni = matRad_calcSigmaIni(this.machine.data,ray,ray.SSD);
            
            % just use tissue classes of voxels found by ray tracer
            if this.calcBioDose
                ray.vTissueIndex = this.vTissueIndex(ray.ix,:);
            end
        end

        function dij = fillDij(this,bixel,dij,stf,currBeamIdx,currRayIdx,currBixelIdx,bixelCounter)
            dij = this.fillDij@DoseEngines.matRad_PencilBeamEngineAbstract(bixel,dij,stf,currBeamIdx,currRayIdx,currBixelIdx,bixelCounter);
            
            % Add MU information
            if ~this.calcDoseDirect
                dij.minMU(bixelCounter,1) = bixel.minMU;
                dij.maxMU(bixelCounter,1) = bixel.maxMU;
                dij.numParticlesPerMU(bixelCounter,1) = bixel.numParticlesPerMU;
            end
        end

        function lateralRayCutOff = getLateralDistanceFromDoseCutOffOnRay(this,ray)
                        % find index of maximum used energy (round to keV for numerical
            % reasons
            maxEnergyIx = max(this.round2(ray.energy,4)) == this.round2([this.machine.data.energy],4);

            lateralRayCutOff = max(this.machine.data(maxEnergyIx).LatCutOff.CutOff);
        end

        function r2 = round2(~,a,b)
            % helper function for energy selection
            r2 = round(a*10^b)/10^b; 
        end

    end   
end

