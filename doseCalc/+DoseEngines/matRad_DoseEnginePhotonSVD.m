classdef matRad_DoseEnginePhotonSVD < DoseEngines.matRad_DoseEnginePencilBeam
% matRad_PhotonDoseEngine: Implements an engine for photon based dose calculation
%   For detailed information see superclass matRad_DoseEngine
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    properties (Constant)
       possibleRadiationModes = 'photons' %constant which represent available radiation modes
       name = 'SVD Pencil Beam'; 
       
       % Define function_Di for beamlet calculation. Constant for use in
       % static computations
       %func_Di = @(x,m,beta) beta/(beta-m) * (exp(-m*x) - exp(-beta*x));
       %func_DiVec = @(x,m,betas) betas./(betas-m) .* (exp(-m*x) - exp(-betas.*x));
    end
    
    properties (SetAccess = public, GetAccess = public)
        useCustomPrimaryPhotonFluence; %boolean to control usage of the primary fluence during dose (influence matrix) computation
        geometricCutOff; %geometrical lateral cutoff [mm]
        kernelCutOff; %cut off in [mm] of kernel values
        randomSeed = 0; %for bixel sampling
        intConvResolution = 0.5; %resolution for kernel convolution [mm]
    end
    
    %Calculation variables
    properties (SetAccess = protected,GetAccess = public)
        isFieldBasedDoseCalc; %Will be set
        penumbraFWHM; %will be obtained from machine
        effectiveLateralCutOff; %will be computed from kernel and geometric lateral cutoffs
        fieldWidth; %Will be obtained during calculation
        
        %Kernel Grid for convolution
        kernelConvSize; %size of the convolution kernel
        kernelX; %meshgrid in X
        kernelZ; %meshgrid in Z
        kernelMxs; %cell array of kernel matrices
        
        gaussFilter; %two-dimensional gaussian filter to model penumbra
        gaussConvSize; %size of the gaussian convolution kernel
        
        convMx_X; %convolution meshgrid in X
        convMx_Z; %convolution meshgrid in Z   
        
        F_X; %fluence meshgrid in X
        F_Z; %fluence meshgrid in Z
        
        
    end
    
   
    methods
        
        function this = matRad_DoseEnginePhotonSVD(ct,stf,pln,cst)
            % Constructor
            % 
            % call
            %   engine = DoseEngines.matRad_DoseEnginePhotonSVD(ct,stf,pln,cst)
            %
            % input
            %   ct:                         matRad ct struct
            %   stf:                        matRad steering information struct
            %   pln:                        matRad plan meta information struct
            %   cst:                        matRad cst struct
                        
            % create this from superclass
            this = this@DoseEngines.matRad_DoseEnginePencilBeam();
            
            if exist('pln','var')
                % 0 if field calc is bixel based, 1 if dose calc is field based
                % num2str is only used to prevent failure of strcmp when bixelWidth
                % contains a number and not a string
                this.isFieldBasedDoseCalc = strcmp(num2str(pln.propStf.bixelWidth),'field');
            end
        end
        
        function dij = calcDose(this,ct,stf,pln,cst)
            % matRad photon dose calculation wrapper
            % can be automaticly called through matRad_calcDose or
            % matRad_calcPhotonDose
            % 
            % call
            %   dij = calcDose(ct,stf,pln,cst)
            %
            % input
            %   ct:             ct cube
            %   stf:            matRad steering information struct
            %   pln:            matRad plan meta information struct
            %   cst:            matRad cst struct
            %   calcDoseDirect: boolian switch to bypass dose influence matrix
            %                   computation and directly calculate dose; only makes
            %                   sense in combination with matRad_calcDoseDirect.m
            %
            % output
            %   dij:            matRad dij struct
            %
            % References
            %   [1] http://www.ncbi.nlm.nih.gov/pubmed/8497215
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


            matRad_cfg =  MatRad_Config.instance();
            matRad_cfg.dispInfo('matRad: Photon dose calculation...\n');
            
            % initialize waitbar           
            figureWait = waitbar(0,'calculate dose influence matrix for photons...');
            
             % show busy state
            set(figureWait,'pointer','watch');
            
            % initialize
            [ct,stf,pln,dij] = this.calcDoseInit(ct,cst,pln,stf);
            
            
            % Precompute kernel convolution if we use a uniform fluence
            if ~this.isFieldBasedDoseCalc   
                % Create fluence matrix
                F = ones(floor(this.fieldWidth/this.intConvResolution));

                if ~this.useCustomPrimaryPhotonFluence
                % gaussian convolution of field to model penumbra
                    F = real(ifft2(fft2(F,this.gaussConvSize,this.gaussConvSize).*fft2(this.gaussFilter,this.gaussConvSize,this.gaussConvSize)));     
                end
            end
           

            counter = 0;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = 1:dij.numOfBeams % loop over all beams
                
                dij = this.calcDoseInitBeam(ct,stf,dij,i);                

                % convolution here if no custom primary fluence and no field based dose calc
                if ~this.useCustomPrimaryPhotonFluence && ~this.isFieldBasedDoseCalc

                    % Display console message.
                    matRad_cfg.dispInfo('\tUniform primary photon fluence -> pre-compute kernel convolution...\n');   
                    
                    % Get kernel interpolators
                    interpKernels = this.getKernelInterpolators(F);                    
                end

                for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray! For field based dose calc, a ray equals a shape

                    counter = counter + 1;
                    this.bixelsPerBeam = this.bixelsPerBeam + 1;

                    % convolution here if custom primary fluence OR field based dose calc
                    if this.useCustomPrimaryPhotonFluence || this.isFieldBasedDoseCalc

                        % overwrite field opening if necessary
                        if this.isFieldBasedDoseCalc
                            F = stf(i).ray(j).shape;
                        end

                        % prepare primary fluence array
                        primaryFluence = this.machine.data.primaryFluence;
                        r     = sqrt( (this.F_X-stf(i).ray(j).rayPos(1)).^2 + (this.F_Z-stf(i).ray(j).rayPos(3)).^2 );
                        Psi   = interp1(primaryFluence(:,1)',primaryFluence(:,2)',r,'linear',0);

                        % apply the primary fluence to the field
                        Fx = F .* Psi;

                        % convolve with the gaussian
                        Fx = real( ifft2(fft2(Fx,this.gaussConvSize,this.gaussConvSize).* fft2(this.gaussFilter,this.gaussConvSize,this.gaussConvSize)) );

                        % Get kernel interpolators
                        interpKernels = this.getKernelInterpolators(Fx);

                    end

                    % Display progress and update text only 200 times
                    if mod(this.bixelsPerBeam,max(1,round(stf(i).totalNumOfBixels/200))) == 0
                        matRad_progress(this.bixelsPerBeam/max(1,round(stf(i).totalNumOfBixels/200)),...
                                        floor(stf(i).totalNumOfBixels/max(1,round(stf(i).totalNumOfBixels/200))));
                    end
                    % update waitbar only 100 times
                    if mod(counter,round(dij.totalNumOfBixels/100)) == 0 && ishandle(figureWait)
                        waitbar(counter/dij.totalNumOfBixels);
                    end

                    % remember beam and bixel number
                    if ~this.calcDoseDirect
                       dij.beamNum(counter)  = i;
                       dij.rayNum(counter)   = j;
                       dij.bixelNum(counter) = 1;
                    else
                        k = 1;
                    end

                    % Ray tracing for beam i and bixel j
                    [ix,rad_distancesSq,isoLatDistsX,isoLatDistsZ] = this.calcGeoDists(this.rot_coordsVdoseGrid, ...
                        stf(i).sourcePoint_bev, ...
                        stf(i).ray(j).targetPoint_bev, ...
                        this.machine.meta.SAD, ...
                        find(~isnan(this.radDepthVdoseGrid{1})), ...
                        this.effectiveLateralCutoff);

                    % empty bixels may happen during recalculation of error
                    % scenarios -> skip to next bixel
                    if isempty(ix)
                        continue;
                    end

                    % calculate photon dose for beam i and bixel j
                    bixelDose = this.calcBixel(interpKernels,ix,isoLatDistsX,isoLatDistsZ);

                    % sample dose only for bixel based dose calculation
                    if ~this.isFieldBasedDoseCalc
                        r0   = 20 + stf(i).bixelWidth;   % [mm] sample beyond the inner core
                        Type = 'radius';
                        [ix,bixelDose] = this.dijSampling(ix,bixelDose,this.radDepthVdoseGrid{1}(ix),rad_distancesSq,Type,r0);
                    end

                    % Save dose for every bixel in cell array
                    this.doseTmpContainer{mod(counter-1,this.numOfBixelsContainer)+1,1} = sparse(this.VdoseGrid(ix),1,bixelDose,dij.doseGrid.numOfVoxels,1);

                    % save computation time and memory 
                    % by sequentially filling the sparse matrix dose.dij from the cell array
                    if mod(counter,this.numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels

                        if this.calcDoseDirect
                            
                            dij = this.fillDijDirect(dij,stf,pln,i,j,k);
                            
                        else
                            
                            dij = this.fillDij(dij,stf,pln,counter);                          

                        end

                    end

                end
                
            end

            %Close Waitbar
            if ishandle(figureWait)
                delete(figureWait);
            end

        end
    end
    
    methods (Access = protected)
        
        function [ct,stf,pln,dij] = calcDoseInit(this,ct,stf,pln,cst)
            %% Assign parameters
            matRad_cfg = MatRad_Config.instance();
            
            % issue warning if biological optimization not possible
            if sum(strcmp(pln.propOpt.bioOptimization,{'effect','RBExD'}))>0
                warndlg('Effect based and RBE optimization not available for photons - physical optimization is carried out instead.');
                pln.bioOptimization = 'none';
            end
            
            % set lateral cutoff value
            if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'geometricCutOff')
                pln.propDoseCalc.geometricCutOff =  matRad_cfg.propDoseCalc.defaultGeometricCutOff; % [mm]
            end

            this.geometricCutOff = pln.propDoseCalc.geometricCutOff;

            if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'kernelCutOff')
                pln.propDoseCalc.kernelCutOff =  matRad_cfg.propDoseCalc.defaultKernelCutOff; % [mm]
            end

            % set kernel cutoff value (determines how much of the kernel is used. This
            % value is separated from lateralCutOff to obtain accurate large open fields)
            this.kernelCutOff = pln.propDoseCalc.kernelCutOff;                      

            % toggle custom primary fluence on/off. if 0 we assume a homogeneous
            % primary fluence, if 1 we use measured radially symmetric data
            if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'useCustomPrimaryPhotonFluence')
                this.useCustomPrimaryPhotonFluence = matRad_cfg.propDoseCalc.defaultUseCustomPrimaryPhotonFluence;
            else
                this.useCustomPrimaryPhotonFluence = pln.propDoseCalc.useCustomPrimaryPhotonFluence;
            end


            % 0 if field calc is bixel based, 1 if dose calc is field based
            % num2str is only used to prevent failure of strcmp when bixelWidth
            % contains a number and not a string
            this.isFieldBasedDoseCalc = strcmp(num2str(pln.propStf.bixelWidth),'field');
            
            %% Call Superclass init
            [ct,stf,pln,dij] = calcDoseInit@DoseEngines.matRad_DoseEnginePencilBeam(this,ct,stf,pln,cst);                        
            
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
            
            if this.kernelCutOff < this.geometricCutOff
                matRad_cfg.dispWarning('Kernel Cut-Off ''%f mm'' cannot be smaller than geometric lateral cutoff ''%f mm''. Using ''%f mm''!',this.kernelCutOff,this.geometricCutOff,this.geometricCutOff);
                this.kernelCutOff = this.geometricCutOff;
            end                        
            
            %% kernel convolution
            % set up convolution grid
            if this.isFieldBasedDoseCalc
                % get data from DICOM import
                this.intConvResolution = pln.propStf.collimation.convResolution; %overwrite default value from dicom
                this.fieldWidth = pln.propStf.collimation.fieldWidth;           
            else
                this.fieldWidth = pln.propStf.bixelWidth;
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
            this.effectiveLateralCutoff = this.geometricCutOff + this.fieldWidth/sqrt(2);
            
            
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
        
        function dij = calcDoseInitBeam(this,ct,stf,dij,i)
            % Method for initializing the beams for analytical pencil beam
            % dose calculation
            %
            % call
            %   this.calcDoseInitBeam(ct,stf,dij,i)
            %
            % input
            %   ct:                         matRad ct struct
            %   stf:                        matRad steering information struct
            %   dij:                        matRad dij struct
            %   i:                          index of beam
            %
            % output
            %   dij:                        updated dij struct
            
            dij = calcDoseInitBeam@DoseEngines.matRad_DoseEnginePencilBeam(this,ct,stf,dij,i);
            
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
        end
            
        
        function dij = fillDij(this,dij,stf,pln,counter)
        % Sequentially fill the sparse matrix dij from the tmpContainer cell arra
        %
        %   see also fillDijDirect
            
            if ~this.calcDoseDirect          
                dij.physicalDose{1}(:,(ceil(counter/this.numOfBixelsContainer)-1)*this.numOfBixelsContainer+1:counter) = [this.doseTmpContainer{1:mod(counter-1,this.numOfBixelsContainer)+1,1}];        
            else
                error([dbstack(1).name ' is not intended for direct dose calculation. For filling the dij inside a direct dose calculation please refer to this.fillDijDirect.']);
            end    
            
        end
        
        function dij = fillDijDirect(this,dij,stf,pln,currBeamIdx,currRayIdx,currBixelIdx)
        % fillDijDirect - sequentially fill dij, meant for direct calculation only
        %   Fill the sparse matrix physicalDose inside dij with the
        %   indices given by the direct dose calculation
        %   
        %   see also fillDij.      
            if this.calcDoseDirect
                if isfield(stf(1).ray(1),'weight') && numel(stf(currBeamIdx).ray(currRayIdx).weight) >= currBixelIdx

                    % score physical dose
                    dij.physicalDose{1}(:,currBeamIdx) = dij.physicalDose{1}(:,currBeamIdx) + stf(currBeamIdx).ray(currRayIdx).weight(currBixelIdx) * this.doseTmpContainer{1,1};
                    
                else
                    error(['No weight available for beam ' num2str(currBeamIdx) ', ray ' num2str(currRayIdx) ', bixel ' num2str(currBixelIdx)]);

                end
            else
                error([dbstack(1).name 'not available for not direct dose calculation. Refer to this.fillDij() for a not direct dose calculation.'])
            end
        end
        
        
        function dose = calcBixel(this,interpKernels,voxelIx,isoLatDistsX,isoLatDistsZ)
            % matRad photon dose calculation for an individual bixel
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
            
            %Here, we just forward to the static implementation
            dose = this.calcSingleBixel(this.machine.meta.SAD,...
                this.machine.data.m,...
                this.machine.data.betas,...
                interpKernels,...
                this.radDepthVdoseGrid{1}(voxelIx),...
                this.geoDistVdoseGrid{1}(voxelIx),...
                isoLatDistsX,...
                isoLatDistsZ);
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
                checkModality = any(strcmp(DoseEngines.matRad_DoseEnginePhotonSVD.possibleRadiationModes, machine.meta.radiationMode));
                
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

