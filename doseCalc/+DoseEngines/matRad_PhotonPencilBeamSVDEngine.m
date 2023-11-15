classdef matRad_PhotonPencilBeamSVDEngine < DoseEngines.matRad_PencilBeamEngineAbstract
    % matRad_PhotonPencilBeamSVDEngine: Pencil-beam dose calculation with 
    % singular value decomposed kernels
    % 
    %
    % References
    %   [1] http://www.ncbi.nlm.nih.gov/pubmed/8497215
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2022 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Constant)
        possibleRadiationModes = {'photons'}; %constant which represent available radiation modes
        name = 'SVD Pencil Beam';
        %supportedQuantities = {'physicalDose'};

        % Define function_Di for beamlet calculation. Constant for use in
        % static computations
        %func_Di = @(x,m,beta) beta/(beta-m) * (exp(-m*x) - exp(-beta*x));
        %func_DiVec = @(x,m,betas) betas./(betas-m) .* (exp(-m*x) - exp(-betas.*x));
    end

    properties (SetAccess = public, GetAccess = public)
        useCustomPrimaryPhotonFluence;  %boolean to control usage of the primary fluence during dose (influence matrix) computation
        kernelCutOff;                   %cut off in [mm] of kernel values
        randomSeed = 0;                 %for bixel sampling
        intConvResolution = 0.5;        %resolution for kernel convolution [mm]

        enableDijSampling = true;
        dijSampling;             %struct with lateral dij sampling parameters
    end

    %Calculation variables
    properties (SetAccess = protected,GetAccess = public)
        isFieldBasedDoseCalc;           %Will be set
        penumbraFWHM;                   %will be obtained from machine
        fieldWidth;                     %Will be obtained during calculation

        %Kernel Grid for convolution
        kernelConvSize;                 %size of the convolution kernel
        kernelX;                        %meshgrid in X
        kernelZ;                        %meshgrid in Z
        kernelMxs;                      %cell array of kernel matrices

        gaussFilter;                    %two-dimensional gaussian filter to model penumbra
        gaussConvSize;                  %size of the gaussian convolution kernel

        convMx_X;                       %convolution meshgrid in X
        convMx_Z;                       %convolution meshgrid in Z

        F_X;                            %fluence meshgrid in X
        F_Z;                            %fluence meshgrid in Z

        Fpre;                           %precomputed fluence if uniform fluence used for calculation
        interpKernelCache;              %Kernel interpolators (cached if precomputation per beam possible)
        
        collimation;                    %collimation structure from dicom import
    end


    methods

        function this = matRad_PhotonPencilBeamSVDEngine(pln)
            % Constructor
            %
            % call
            %   engine = DoseEngines.matRad_PhotonPencilBeamSVDEngine(ct,stf,pln,cst)
            %
            % input
            %   ct:                         matRad ct struct
            %   stf:                        matRad steering information struct
            %   pln:                        matRad plan meta information struct
            %   cst:                        matRad cst struct

            % create this from superclass
            this = this@DoseEngines.matRad_PencilBeamEngineAbstract(pln);            

            if nargin > 0
                % 0 if field calc is bixel based, 1 if dose calc is field based
                % num2str is only used to prevent failure of strcmp when bixelWidth
                % contains a number and not a string
                this.isFieldBasedDoseCalc = strcmp(num2str(pln.propStf.bixelWidth),'field');
            end
        end

        function setDefaults(this)
            setDefaults@DoseEngines.matRad_PencilBeamEngineAbstract(this);

            %Assign defaults from Config
            matRad_cfg = MatRad_Config.instance();
            this.useCustomPrimaryPhotonFluence  = matRad_cfg.propDoseCalc.defaultUseCustomPrimaryPhotonFluence;
            this.kernelCutOff                   = matRad_cfg.propDoseCalc.defaultKernelCutOff;
            
            %dij sampling defaults                      
            this.dijSampling.relDoseThreshold = 0.01;
            this.dijSampling.latCutOff        = 20; 
            this.dijSampling.type             = 'radius';
            this.dijSampling.deltaRadDepth    = 5;
        end
    end

    methods (Access = protected)

        function [dij,ct,cst,stf] = calcDoseInit(this,ct,cst,stf)
            %% Assign parameters
            matRad_cfg = MatRad_Config.instance();

            % 0 if field calc is bixel based, 1 if dose calc is field based
            % num2str is only used to prevent failure of strcmp when bixelWidth
            % contains a number and not a string
            this.isFieldBasedDoseCalc = any(arrayfun(@(s) strcmp(num2str(s.bixelWidth),'field'),stf));

            %% Call Superclass init
            [dij,ct,cst,stf] = calcDoseInit@DoseEngines.matRad_PencilBeamEngineAbstract(this,ct,cst,stf);

            %% Validate some properties
            % gaussian filter to model penumbra from (measured) machine output / see
            % diploma thesis siggel 4.1.2 -> https://github.com/e0404/matRad/wiki/Dose-influence-matrix-calculation
            if isfield(this.machine.data,'penumbraFWHMatIso')
                this.penumbraFWHM = this.machine.data.penumbraFWHMatIso;
            else
                this.penumbraFWHM = 5;
                matRad_cfg.dispWarning('photon machine file does not contain measured penumbra width in machine.data.penumbraFWHMatIso. Assuming %f mm.',this.penumbraFWHM);
            end

            %Correct kernel cut off to base data limits if needed
            if this.kernelCutOff > this.machine.data.kernelPos(end)
                matRad_cfg.dispWarning('Kernel Cut-Off ''%f mm'' larger than machine data range of ''%f mm''. Using ''%f mm''!',this.kernelCutOff,this.machine.data.kernelPos(end),this.machine.data.kernelPos(end));
                this.kernelCutOff = this.machine.data.kernelPos(end);
            end

            if this.kernelCutOff < this.geometricLateralCutOff
                matRad_cfg.dispWarning('Kernel Cut-Off ''%f mm'' cannot be smaller than geometric lateral cutoff ''%f mm''. Using ''%f mm''!',this.kernelCutOff,this.geometricLateralCutOff,this.geometricLateralCutOff);
                this.kernelCutOff = this.geometricLateralCutOff;
            end

            %% kernel convolution
            % set up convolution grid
            if this.isFieldBasedDoseCalc
                % get data from DICOM import
                this.intConvResolution = this.collimation.convResolution; %overwrite default value from dicom
                this.fieldWidth = this.collimation.fieldWidth;
            else
                if numel(unique([stf.bixelWidth])) > 1
                    matRad_cfg.dispError('Different bixelWidths pear beam are not supported!');
                end

                this.fieldWidth = unique([stf.bixelWidth]);
            end

            % calculate field size and distances
            fieldLimit = ceil(this.fieldWidth/(2*this.intConvResolution));
            [this.F_X,this.F_Z] = meshgrid(-fieldLimit*this.intConvResolution: ...
                this.intConvResolution: ...
                (fieldLimit-1)*this.intConvResolution);



            sigmaGauss = this.penumbraFWHM / sqrt(8*log(2)); % [mm]
            % use 5 times sigma as the limits for the gaussian convolution
            gaussLimit = ceil(5*sigmaGauss/this.intConvResolution);
            [gaussFilterX,gaussFilterZ] = meshgrid(-gaussLimit*this.intConvResolution: ...
                this.intConvResolution: ...
                (gaussLimit-1)*this.intConvResolution);
            this.gaussFilter =  1/(2*pi*sigmaGauss^2/this.intConvResolution^2) * exp(-(gaussFilterX.^2+gaussFilterZ.^2)/(2*sigmaGauss^2) );
            this.gaussConvSize = 2*(fieldLimit + gaussLimit);

            % get kernel size and distances

            kernelLimit = ceil(this.kernelCutOff/this.intConvResolution);
            [this.kernelX, this.kernelZ] = meshgrid(-kernelLimit*this.intConvResolution: ...
                this.intConvResolution: ...
                (kernelLimit-1)*this.intConvResolution);

            % precalculate convolved kernel size and distances
            kernelConvLimit = fieldLimit + gaussLimit + kernelLimit;
            [this.convMx_X, this.convMx_Z] = meshgrid(-kernelConvLimit*this.intConvResolution: ...
                this.intConvResolution: ...
                (kernelConvLimit-1)*this.intConvResolution);
            % calculate also the total size and distance as we need this during convolution extensively
            this.kernelConvSize = 2*kernelConvLimit;

            % define an effective lateral cutoff where dose will be calculated. note
            % that storage within the influence matrix may be subject to sampling
            this.effectiveLateralCutOff = this.geometricLateralCutOff + this.fieldWidth/sqrt(2);
            

            % Check if we can precompute fluence and precompute kernel 
            % convolution if we use a uniform fluence
            if ~this.isFieldBasedDoseCalc
                % Create fluence matrix
                this.Fpre = ones(floor(this.fieldWidth/this.intConvResolution));

                if ~this.useCustomPrimaryPhotonFluence
                    % gaussian convolution of field to model penumbra
                    this.Fpre = real(ifft2(fft2(this.Fpre,this.gaussConvSize,this.gaussConvSize).*fft2(this.gaussFilter,this.gaussConvSize,this.gaussConvSize)));
                end
            end

            %% Initialize randomization
            [env, ~] = matRad_getEnvironment();

            switch env
                case 'MATLAB'
                    rng(this.randomSeed); %Initializes Mersenne Twister with seed 0
                case 'OCTAVE'
                    rand('state',this.randomSeed); %Initializes Mersenne Twister with state 0 (does not give similar random numbers as in Matlab)
                otherwise
                    rand('seed',this.randomSeed); %Fallback
                    matRad_cfg.dispWarning('Environment %s not recognized!',env);
            end

        end

        function currBeam = initBeam(this,currBeam,ct,cst,stf,i)
            % Method for initializing the beams for analytical pencil beam
            % dose calculation
            %
            % call
            %   this.initBeam(ct,stf,dij,i)
            %
            % input
            %   ct:                         matRad ct struct
            %   stf:                        matRad steering information struct
            %   dij:                        matRad dij struct
            %   i:                          index of beam
            %
            % output
            %   dij:                        updated dij struct

            currBeam = initBeam@DoseEngines.matRad_PencilBeamEngineAbstract(this,currBeam,ct,cst,stf,i);
            
            currBeam.ixRadDepths = find( ~isnan(currBeam.radDepthVdoseGrid{1}));
            currBeam.rot_coordsVdoseGrid = currBeam.rot_coordsVdoseGrid(currBeam.ixRadDepths,:);

            matRad_cfg = MatRad_Config.instance();

            % get index of central ray or closest to the central ray
            [~,center] = min(sum(reshape([stf(i).ray.rayPos_bev],3,[]).^2));

            % get correct kernel for given SSD at central ray (nearest neighbor approximation)
            [~,currSSDix] = min(abs([this.machine.data.kernel.SSD]-stf(i).ray(center).SSD));
            % Display console message.
            matRad_cfg.dispInfo('\tSSD = %g mm ...\n',this.machine.data.kernel(currSSDix).SSD);

            %Hardcoded for now
            useKernels = {'kernel1','kernel2','kernel3'};

            kernelPos = this.machine.data.kernelPos;

            for k = 1:length(useKernels)
                kernel = this.machine.data.kernel(currSSDix).(useKernels{k});
                this.kernelMxs{k} = interp1(kernelPos,kernel,sqrt(this.kernelX.^2+this.kernelZ.^2),'linear',0);
            end

            % convolution here if no custom primary fluence and no field based dose calc
            if ~isempty(this.Fpre) && ~this.useCustomPrimaryPhotonFluence && ~this.isFieldBasedDoseCalc

                % Display console message.
                matRad_cfg.dispInfo('\tUniform primary photon fluence -> pre-compute kernel convolution...\n');

                % Get kernel interpolators
                this.interpKernelCache = this.getKernelInterpolators(this.Fpre);
            end
        end

        function [bixel] = computeBixel(this,currRay,k)
            % matRad photon dose calculation for an individual bixel
            %
            % call
            %   bixel = this.computeBixel(currRay,k)
            
            bixel = struct();

            if isfield(this.tmpMatrixContainers,'physicalDose')            
                bixel.physicalDose = this.calcSingleBixel(currRay.SAD,...
                    this.machine.data.m,...
                    this.machine.data.betas,...
                    currRay.interpKernels,...
                    currRay.radDepths,...
                    currRay.geoDepths,...
                    currRay.isoLatDists(:,1),...
                    currRay.isoLatDists(:,2));

                % sample dose only for bixel based dose calculation
                if this.enableDijSampling && ~this.isFieldBasedDoseCalc
                    [bixel.ix,bixel.physicalDose] = this.sampleDij(currRay.ix,bixel.physicalDose,currRay.radDepths,currRay.radialDist_sq,currRay.bixelWidth);
                else
                    bixel.ix = currRay.ix;
                end
            else
                bixel.ix = [];
            end
        end
        
        function interpKernels = getKernelInterpolators(this,Fx)

            matRad_cfg = MatRad_Config.instance();

            nKernels = length(this.kernelMxs);
            interpKernels = cell(1,nKernels);

            for ik = 1:nKernels
                % 2D convolution of Fluence and Kernels in fourier domain
                convMx = real( ifft2(fft2(Fx,this.kernelConvSize,this.kernelConvSize).* fft2(this.kernelMxs{ik},this.kernelConvSize,this.kernelConvSize)));

                % Creates an interpolant for kernes from vectors position X and Z
                if matRad_cfg.isMatlab
                    interpKernels{ik} = griddedInterpolant(this.convMx_X',this.convMx_Z',convMx','linear','none');
                elseif matRad_cfg.isOctave
                    %For some reason the use of interpn here is much faster
                    %than using interp2 in Octave
                    interpKernels{ik} = @(x,y) interpn(this.convMx_X(1,:),this.convMx_Z(:,1),convMx',x,y,'linear',NaN);
                end
            end
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
                LatCutOff = LatCutOff + bixelWidth/sqrt(2); %use half of the bixel width diagonal as max. field size radius for sampling
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
    
        function [ray] = initRay(this,currBeam,j)

            ray = initRay@DoseEngines.matRad_PencilBeamEngineAbstract(this,currBeam,j);

            % convolution here if custom primary fluence OR field based dose calc
            if this.useCustomPrimaryPhotonFluence || this.isFieldBasedDoseCalc

                % overwrite field opening if necessary
                if this.isFieldBasedDoseCalc
                    F = ray.shape;
                end

                % prepare primary fluence array
                primaryFluence = this.machine.data.primaryFluence;
                r     = sqrt( (this.F_X-ray.rayPos(1)).^2 + (this.F_Z-ray.rayPos(3)).^2 );
                Psi   = interp1(primaryFluence(:,1)',primaryFluence(:,2)',r,'linear',0);

                % apply the primary fluence to the field
                Fx = F .* Psi;

                % convolve with the gaussian
                Fx = real( ifft2(fft2(Fx,this.gaussConvSize,this.gaussConvSize).* fft2(this.gaussFilter,this.gaussConvSize,this.gaussConvSize)) );

                % Get kernel interpolators
                ray.interpKernels = this.getKernelInterpolators(Fx);

            else
                ray.interpKernels = this.interpKernelCache;
            end

            ray.geoDepths = currBeam.geoDistVdoseGrid{1}(ray.ix);
        end    
    end

    methods (Static)

        function [available,msg] = isAvailable(pln,machine)
            % see superclass for information

            msg = [];
            available = false;

            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            %checkBasic
            try
                checkBasic = isfield(machine,'meta') && isfield(machine,'data');

                %check modality
                checkModality = any(strcmp(DoseEngines.matRad_PhotonPencilBeamSVDEngine.possibleRadiationModes, machine.meta.radiationMode));

                preCheck = checkBasic && checkModality;

                if ~preCheck
                    return;
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return;
            end


            %Basic check for information (does not check data integrity & subfields etc.)
            checkData = all(isfield(machine.data,{'betas','energy','m','primaryFluence','kernel','kernelPos'}));
            checkMeta = all(isfield(machine.meta,{'SAD','SCD'}));

            if checkData && checkMeta
                available = true;
            else
                available = false;
                return;
            end

            %Now check for optional fields that would be guessed otherwise
            checkOptional = isfield(machine.data,'penumbraFWHMatIso');
            if checkOptional
                msg = 'No penumbra given, generic value will be used!';
            end
        end

        function bixelDose = calcSingleBixel(SAD,m,betas,interpKernels,...
                radDepths,geoDists,isoLatDistsX,isoLatDistsZ)
            % matRad photon dose calculation for an individual bixel
            %   This is defined as a static function so it can also be
            %   called individually for certain applications without having
            %   a fully defined dose engine
            %
            % call
            %   dose = this.calcPhotonDoseBixel(SAD,m,betas,Interp_kernel1,...
            %                  Interp_kernel2,Interp_kernel3,radDepths,geoDists,...
            %                  isoLatDistsX,isoLatDistsZ)
            %
            % input
            %   SAD:                source to axis distance
            %   m:                  absorption in water (part of the dose calc base
            %                       data)
            %   betas:              beta parameters for the parameterization of the
            %                       three depth dose components
            %   interpKernels:      kernel interpolators for dose calculation
            %   radDepths:          radiological depths
            %   geoDists:           geometrical distance from virtual photon source
            %   isoLatDistsX:       lateral distance in X direction in BEV from central
            %                       ray at iso center plane
            %   isoLatDistsZ:       lateral distance in Z direction in BEV from central
            %                       ray at iso center plane
            %
            % output
            %   dose:               photon dose at specified locations as linear vector
            %
            % References
            %   [1] http://www.ncbi.nlm.nih.gov/pubmed/8497215
            %

            % Compute depth dose components according to [1, eq. 17]
            doseComponent = betas./(betas-m) .* (exp(-m*radDepths) - exp(-betas.*radDepths));

            % Multiply with lateral 2D-convolved kernels using
            % grid interpolation at lateral distances (summands in [1, eq.
            % 19] w/o inv sq corr)
            for ik = 1:length(interpKernels)
                doseComponent(:,ik) = doseComponent(:,ik) .* interpKernels{ik}(isoLatDistsX,isoLatDistsZ);
            end

            % now add everything together (eq 19 w/o inv sq corr -> see below)
            bixelDose = sum(doseComponent,2);

            % inverse square correction
            bixelDose = bixelDose .* ((SAD)./geoDists(:)).^2;

            % check if we have valid dose values and adjust numerical instabilities
            % from fft convolution
            bixelDose(bixelDose < 0 & bixelDose > -1e-14) = 0;
            if any(isnan(bixelDose)) || any(bixelDose<0)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid numerical values in photon dose calculation.');
            end
        end

    end

end

