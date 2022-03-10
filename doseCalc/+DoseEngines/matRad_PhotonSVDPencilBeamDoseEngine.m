classdef matRad_PhotonSVDPencilBeamDoseEngine < DoseEngines.matRad_AnalyticalPencilBeamEngine
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
        
        function obj = matRad_PhotonSVDPencilBeamDoseEngine(ct,stf,pln,cst)
            % Constructor
            % 
            % call
            %   engine = DoseEngines.matRad_PhotonSVDPencilBeamDoseEngine(ct,stf,pln,cst)
            %
            % input
            %   ct:                         matRad ct struct
            %   stf:                        matRad steering information struct
            %   pln:                        matRad plan meta information struct
            %   cst:                        matRad cst struct
                        
            % create obj of superclass
            obj = obj@DoseEngines.matRad_AnalyticalPencilBeamEngine();
            
            if exist('pln','var')
                % 0 if field calc is bixel based, 1 if dose calc is field based
                % num2str is only used to prevent failure of strcmp when bixelWidth
                % contains a number and not a string
                obj.isFieldBasedDoseCalc = strcmp(num2str(pln.propStf.bixelWidth),'field');
            end
        end
        
        function dij = calcDose(obj,ct,stf,pln,cst)
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
            [ct,stf,pln,dij] = obj.calcDoseInit(ct,stf,pln,cst);
            
            
            % Precompute kernel convolution if we use a uniform fluence
            if ~obj.isFieldBasedDoseCalc   
                % Create fluence matrix
                F = ones(floor(obj.fieldWidth/obj.intConvResolution));

                if ~obj.useCustomPrimaryPhotonFluence
                % gaussian convolution of field to model penumbra
                    F = real(ifft2(fft2(F,obj.gaussConvSize,obj.gaussConvSize).*fft2(obj.gaussFilter,obj.gaussConvSize,obj.gaussConvSize)));     
                end
            end
           

            counter = 0;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = 1:dij.numOfBeams % loop over all beams
                
                dij = obj.calcDoseInitBeam(ct,stf,dij,i);                

                % convolution here if no custom primary fluence and no field based dose calc
                if ~obj.useCustomPrimaryPhotonFluence && ~obj.isFieldBasedDoseCalc

                    % Display console message.
                    matRad_cfg.dispInfo('\tUniform primary photon fluence -> pre-compute kernel convolution...\n');   
                    
                    % Get kernel interpolators
                    interpKernels = obj.getKernelInterpolators(F);                    
                end

                for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray! For field based dose calc, a ray equals a shape

                    counter = counter + 1;
                    obj.bixelsPerBeam = obj.bixelsPerBeam + 1;

                    % convolution here if custom primary fluence OR field based dose calc
                    if obj.useCustomPrimaryPhotonFluence || obj.isFieldBasedDoseCalc

                        % overwrite field opening if necessary
                        if obj.isFieldBasedDoseCalc
                            F = stf(i).ray(j).shape;
                        end

                        % prepare primary fluence array
                        primaryFluence = obj.machine.data.primaryFluence;
                        r     = sqrt( (obj.F_X-stf(i).ray(j).rayPos(1)).^2 + (obj.F_Z-stf(i).ray(j).rayPos(3)).^2 );
                        Psi   = interp1(primaryFluence(:,1)',primaryFluence(:,2)',r,'linear',0);

                        % apply the primary fluence to the field
                        Fx = F .* Psi;

                        % convolve with the gaussian
                        Fx = real( ifft2(fft2(Fx,obj.gaussConvSize,obj.gaussConvSize).* fft2(obj.gaussFilter,obj.gaussConvSize,obj.gaussConvSize)) );

                        % Get kernel interpolators
                        interpKernels = obj.getKernelInterpolators(Fx);

                    end

                    % Display progress and update text only 200 times
                    if mod(obj.bixelsPerBeam,max(1,round(stf(i).totalNumOfBixels/200))) == 0
                        matRad_progress(obj.bixelsPerBeam/max(1,round(stf(i).totalNumOfBixels/200)),...
                                        floor(stf(i).totalNumOfBixels/max(1,round(stf(i).totalNumOfBixels/200))));
                    end
                    % update waitbar only 100 times
                    if mod(counter,round(dij.totalNumOfBixels/100)) == 0 && ishandle(figureWait)
                        waitbar(counter/dij.totalNumOfBixels);
                    end

                    % remember beam and bixel number
                    if ~obj.calcDoseDirect
                       dij.beamNum(counter)  = i;
                       dij.rayNum(counter)   = j;
                       dij.bixelNum(counter) = 1;
                    else
                        k = 1;
                    end

                    % Ray tracing for beam i and bixel j
                    [ix,rad_distancesSq,isoLatDistsX,isoLatDistsZ] = obj.calcGeoDists(obj.rot_coordsVdoseGrid, ...
                        stf(i).sourcePoint_bev, ...
                        stf(i).ray(j).targetPoint_bev, ...
                        obj.machine.meta.SAD, ...
                        find(~isnan(obj.radDepthVdoseGrid{1})), ...
                        obj.effectiveLateralCutoff);

                    % empty bixels may happen during recalculation of error
                    % scenarios -> skip to next bixel
                    if isempty(ix)
                        continue;
                    end

                    % calculate photon dose for beam i and bixel j
                    bixelDose = obj.calcBixel(interpKernels,ix,isoLatDistsX,isoLatDistsZ);

                    % sample dose only for bixel based dose calculation
                    if ~obj.isFieldBasedDoseCalc
                        r0   = 20 + stf(i).bixelWidth;   % [mm] sample beyond the inner core
                        Type = 'radius';
                        [ix,bixelDose] = obj.dijSampling(ix,bixelDose,obj.radDepthVdoseGrid{1}(ix),rad_distancesSq,Type,r0);
                    end

                    % Save dose for every bixel in cell array
                    obj.doseTmpContainer{mod(counter-1,obj.numOfBixelsContainer)+1,1} = sparse(obj.VdoseGrid(ix),1,bixelDose,dij.doseGrid.numOfVoxels,1);

                    % save computation time and memory 
                    % by sequentially filling the sparse matrix dose.dij from the cell array
                    if mod(counter,obj.numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels

                        if obj.calcDoseDirect
                            
                            dij = obj.fillDijDirect(dij,stf,pln,i,j,k);
                            
                        else
                            
                            dij = obj.fillDij(dij,stf,pln,counter);                          

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
        
        function [ct,stf,pln,dij] = calcDoseInit(obj,ct,stf,pln,cst)
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

            obj.geometricCutOff = pln.propDoseCalc.geometricCutOff;

            if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'kernelCutOff')
                pln.propDoseCalc.kernelCutOff =  matRad_cfg.propDoseCalc.defaultKernelCutOff; % [mm]
            end

            % set kernel cutoff value (determines how much of the kernel is used. This
            % value is separated from lateralCutOff to obtain accurate large open fields)
            obj.kernelCutOff = pln.propDoseCalc.kernelCutOff;                      

            % toggle custom primary fluence on/off. if 0 we assume a homogeneous
            % primary fluence, if 1 we use measured radially symmetric data
            if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'useCustomPrimaryPhotonFluence')
                obj.useCustomPrimaryPhotonFluence = matRad_cfg.propDoseCalc.defaultUseCustomPrimaryPhotonFluence;
            else
                obj.useCustomPrimaryPhotonFluence = pln.propDoseCalc.useCustomPrimaryPhotonFluence;
            end


            % 0 if field calc is bixel based, 1 if dose calc is field based
            % num2str is only used to prevent failure of strcmp when bixelWidth
            % contains a number and not a string
            obj.isFieldBasedDoseCalc = strcmp(num2str(pln.propStf.bixelWidth),'field');
            
            %% Call Superclass init
            [ct,stf,pln,dij] = calcDoseInit@DoseEngines.matRad_AnalyticalPencilBeamEngine(obj,ct,stf,pln,cst);                        
            
            %% Validate some properties
            % gaussian filter to model penumbra from (measured) machine output / see
            % diploma thesis siggel 4.1.2 -> https://github.com/e0404/matRad/wiki/Dose-influence-matrix-calculation
            if isfield(obj.machine.data,'penumbraFWHMatIso')
                obj.penumbraFWHM = obj.machine.data.penumbraFWHMatIso;
            else
                obj.penumbraFWHM = 5;
                matRad_cfg.dispWarning('photon machine file does not contain measured penumbra width in machine.data.penumbraFWHMatIso. Assuming %f mm.',obj.penumbraFWHM);
            end
            
            %Correct kernel cut off to base data limits if needed
            if obj.kernelCutOff > obj.machine.data.kernelPos(end)
                matRad_cfg.dispWarning('Kernel Cut-Off ''%f mm'' larger than machine data range of ''%f mm''. Using ''%f mm''!',obj.kernelCutOff,obj.machine.data.kernelPos(end),obj.machine.data.kernelPos(end));
                obj.kernelCutOff = obj.machine.data.kernelPos(end);
            end
            
            if obj.kernelCutOff < obj.geometricCutOff
                matRad_cfg.dispWarning('Kernel Cut-Off ''%f mm'' cannot be smaller than geometric lateral cutoff ''%f mm''. Using ''%f mm''!',obj.kernelCutOff,obj.geometricCutOff,obj.geometricCutOff);
                obj.kernelCutOff = obj.geometricCutOff;
            end                        
            
            %% kernel convolution
            % set up convolution grid
            if obj.isFieldBasedDoseCalc
                % get data from DICOM import
                obj.intConvResolution = pln.propStf.collimation.convResolution; %overwrite default value from dicom
                obj.fieldWidth = pln.propStf.collimation.fieldWidth;           
            else
                obj.fieldWidth = pln.propStf.bixelWidth;
            end

            % calculate field size and distances
            fieldLimit = ceil(obj.fieldWidth/(2*obj.intConvResolution));
            [obj.F_X,obj.F_Z] = meshgrid(-fieldLimit*obj.intConvResolution: ...
                              obj.intConvResolution: ...
                              (fieldLimit-1)*obj.intConvResolution);    

            

            sigmaGauss = obj.penumbraFWHM / sqrt(8*log(2)); % [mm] 
            % use 5 times sigma as the limits for the gaussian convolution
            gaussLimit = ceil(5*sigmaGauss/obj.intConvResolution);
            [gaussFilterX,gaussFilterZ] = meshgrid(-gaussLimit*obj.intConvResolution: ...
                                                obj.intConvResolution: ...
                                                (gaussLimit-1)*obj.intConvResolution);   
            obj.gaussFilter =  1/(2*pi*sigmaGauss^2/obj.intConvResolution^2) * exp(-(gaussFilterX.^2+gaussFilterZ.^2)/(2*sigmaGauss^2) );
            obj.gaussConvSize = 2*(fieldLimit + gaussLimit);

            % get kernel size and distances

            kernelLimit = ceil(obj.kernelCutOff/obj.intConvResolution);     
            [obj.kernelX, obj.kernelZ] = meshgrid(-kernelLimit*obj.intConvResolution: ...
                                        obj.intConvResolution: ...
                                        (kernelLimit-1)*obj.intConvResolution);

            % precalculate convolved kernel size and distances
            kernelConvLimit = fieldLimit + gaussLimit + kernelLimit;
            [obj.convMx_X, obj.convMx_Z] = meshgrid(-kernelConvLimit*obj.intConvResolution: ...
                                            obj.intConvResolution: ...
                                            (kernelConvLimit-1)*obj.intConvResolution);
            % calculate also the total size and distance as we need this during convolution extensively
            obj.kernelConvSize = 2*kernelConvLimit;

            % define an effective lateral cutoff where dose will be calculated. note
            % that storage within the influence matrix may be subject to sampling
            obj.effectiveLateralCutoff = obj.geometricCutOff + obj.fieldWidth/sqrt(2);
            
            
            %% Initialize randomization
            [env, ~] = matRad_getEnvironment();
            
            switch env
                case 'MATLAB'
                    rng(obj.randomSeed); %Initializes Mersenne Twister with seed 0
                case 'OCTAVE'
                    rand('state',obj.randomSeed); %Initializes Mersenne Twister with state 0 (does not give similar random numbers as in Matlab)
                otherwise
                    rand('seed',obj.randomSeed); %Fallback
                    matRad_cfg.dispWarning('Environment %s not recognized!',env);
            end

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
            
            dij = calcDoseInitBeam@DoseEngines.matRad_AnalyticalPencilBeamEngine(obj,ct,stf,dij,i);
            
            matRad_cfg = MatRad_Config.instance();
            
            % get index of central ray or closest to the central ray
            [~,center] = min(sum(reshape([stf(i).ray.rayPos_bev],3,[]).^2));
            
            % get correct kernel for given SSD at central ray (nearest neighbor approximation)
            [~,currSSDix] = min(abs([obj.machine.data.kernel.SSD]-stf(i).ray(center).SSD));
            % Display console message.
            matRad_cfg.dispInfo('\tSSD = %g mm ...\n',obj.machine.data.kernel(currSSDix).SSD);
            
            %Hardcoded for now
            useKernels = {'kernel1','kernel2','kernel3'};
            
            kernelPos = obj.machine.data.kernelPos;
            
            for k = 1:length(useKernels)
                kernel = obj.machine.data.kernel(currSSDix).(useKernels{k});
                obj.kernelMxs{k} = interp1(kernelPos,kernel,sqrt(obj.kernelX.^2+obj.kernelZ.^2),'linear',0);
            end
        end
            
        
        function dij = fillDij(obj,dij,stf,pln,counter)
        % Sequentially fill the sparse matrix dij from the tmpContainer cell arra
        %
        %   see also fillDijDirect
            
            if ~obj.calcDoseDirect          
                dij.physicalDose{1}(:,(ceil(counter/obj.numOfBixelsContainer)-1)*obj.numOfBixelsContainer+1:counter) = [obj.doseTmpContainer{1:mod(counter-1,obj.numOfBixelsContainer)+1,1}];        
            else
                error([dbstack(1).name ' is not intended for direct dose calculation. For filling the dij inside a direct dose calculation please refer to obj.fillDijDirect.']);
            end    
            
        end
        
        function dij = fillDijDirect(obj,dij,stf,pln,currBeamIdx,currRayIdx,currBixelIdx)
        % fillDijDirect - sequentially fill dij, meant for direct calculation only
        %   Fill the sparse matrix physicalDose inside dij with the
        %   indices given by the direct dose calculation
        %   
        %   see also fillDij.      
            if obj.calcDoseDirect
                if isfield(stf(1).ray(1),'weight') && numel(stf(currBeamIdx).ray(currRayIdx).weight) >= currBixelIdx

                    % score physical dose
                    dij.physicalDose{1}(:,currBeamIdx) = dij.physicalDose{1}(:,currBeamIdx) + stf(currBeamIdx).ray(currRayIdx).weight(currBixelIdx) * obj.doseTmpContainer{1,1};
                    
                else
                    error(['No weight available for beam ' num2str(currBeamIdx) ', ray ' num2str(currRayIdx) ', bixel ' num2str(currBixelIdx)]);

                end
            else
                error([dbstack(1).name 'not available for not direct dose calculation. Refer to obj.fillDij() for a not direct dose calculation.'])
            end
        end
        
        
        function dose = calcBixel(obj,interpKernels,voxelIx,isoLatDistsX,isoLatDistsZ)
            % matRad photon dose calculation for an individual bixel
            % 
            % call
            %   dose = obj.calcPhotonDoseBixel(SAD,m,betas,Interp_kernel1,...
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
            dose = obj.calcSingleBixel(obj.machine.meta.SAD,...
                obj.machine.data.m,...
                obj.machine.data.betas,...
                interpKernels,...
                obj.radDepthVdoseGrid{1}(voxelIx),...
                obj.geoDistVdoseGrid{1}(voxelIx),...
                isoLatDistsX,...
                isoLatDistsZ);
        end
        
        
        function interpKernels = getKernelInterpolators(obj,Fx)
            
            matRad_cfg = MatRad_Config.instance();
            
            nKernels = length(obj.kernelMxs);
            interpKernels = cell(1,nKernels);
            
            for ik = 1:nKernels                
                % 2D convolution of Fluence and Kernels in fourier domain
                convMx = real( ifft2(fft2(Fx,obj.kernelConvSize,obj.kernelConvSize).* fft2(obj.kernelMxs{ik},obj.kernelConvSize,obj.kernelConvSize)));
                
                % Creates an interpolant for kernes from vectors position X and Z
                if matRad_cfg.isMatlab
                    interpKernels{ik} = griddedInterpolant(obj.convMx_X',obj.convMx_Z',convMx','linear','none');
                elseif matRad_cfg.isOctave
                    %For some reason the use of interpn here is much faster
                    %than using interp2 in Octave
                    interpKernels{ik} = @(x,y) interpn(obj.convMx_X(1,:),obj.convMx_Z(:,1),convMx',x,y,'linear',NaN);
                end
            end
        end
    end
    
    methods (Static)

        function ret = isAvailable(pln)
            % see superclass for information
            ret = any(strcmp(DoseEngines.matRad_PhotonSVDPencilBeamDoseEngine.possibleRadiationModes, pln.radiationMode));
        end
        
        function bixelDose = calcSingleBixel(SAD,m,betas,interpKernels,...
                radDepths,geoDists,isoLatDistsX,isoLatDistsZ)
            % matRad photon dose calculation for an individual bixel
            %   This is defined as a static function so it can also be
            %   called individually for certain applications without having
            %   a fully defined dose engine
            % 
            % call
            %   dose = obj.calcPhotonDoseBixel(SAD,m,betas,Interp_kernel1,...
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

