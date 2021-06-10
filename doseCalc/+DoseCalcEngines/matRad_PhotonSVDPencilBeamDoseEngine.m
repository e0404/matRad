classdef matRad_PhotonSVDPencilBeamDoseEngine < DoseCalcEngines.matRad_AnalyticalPencilBeamEngine
% matRad_PhotonDoseEngine: Implements an engine for photon based dose calculation
%   For detailed information see superclass matRad_DoseCalcEngine
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
       possibleRadiationModes = "photons" %constant which represent available radiation modes
       name = "photon dose calculation";
    end
    
    properties (SetAccess = private, GetAccess = public)
        isFieldBasedDoseCalc;
    end
    
   
    methods
        
        function obj = matRad_PhotonSVDPencilBeamDoseEngine(ct,stf,pln,cst,calcDoseDirect)
            
            if nargin == 0 || ~exist('calcDoseDirect','var')
                calcDoseDirect = false;
            else
                % future code here 
            end
            
            % create obj of superclass
            obj = obj@DoseCalcEngines.matRad_AnalyticalPencilBeamEngine(calcDoseDirect);
            
            % 0 if field calc is bixel based, 1 if dose calc is field based
            % num2str is only used to prevent failure of strcmp when bixelWidth
            % contains a number and not a string
            obj.isFieldBasedDoseCalc = strcmp(num2str(pln.propStf.bixelWidth),'field');
        end
        
        function dij = calculateDose(obj,ct,stf,pln,cst)
            % matRad photon dose calculation wrapper
            % 
            % call
            %   dij = matRad_calcPhotonDose(ct,stf,pln,cst)
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

            % initialize
            [ct,stf,pln,dij] = obj.calcDoseInit(ct,stf,pln,cst);

            [env, ~] = matRad_getEnvironment();

            switch env
                 case 'MATLAB'
                      rng('default');
                 case 'OCTAVE'
                      rand('seed',0)
            end

            % issue warning if biological optimization not possible
            if sum(strcmp(pln.propOpt.bioOptimization,{'effect','RBExD'}))>0
                warndlg('Effect based and RBE optimization not available for photons - physical optimization is carried out instead.');
                pln.bioOptimization = 'none';
            end

            % initialize waitbar
            figureWait = waitbar(0,'calculate dose influence matrix for photons...');
            % show busy state
            set(figureWait,'pointer','watch');

            % set lateral cutoff value
            lateralCutoff = matRad_cfg.propDoseCalc.defaultGeometricCutOff; % [mm]

            % toggle custom primary fluence on/off. if 0 we assume a homogeneous
            % primary fluence, if 1 we use measured radially symmetric data
            if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'useCustomPrimaryPhotonFluence')
                useCustomPrimFluenceBool = matRad_cfg.propDoseCalc.defaultUseCustomPrimaryPhotonFluence;
            else
                useCustomPrimFluenceBool = pln.propDoseCalc.useCustomPrimaryPhotonFluence;
            end


            % 0 if field calc is bixel based, 1 if dose calc is field based
            % num2str is only used to prevent failure of strcmp when bixelWidth
            % contains a number and not a string
            obj.isFieldBasedDoseCalc = strcmp(num2str(pln.propStf.bixelWidth),'field');

            % kernel convolution
            % set up convolution grid
            if obj.isFieldBasedDoseCalc
                % get data from DICOM import
                intConvResolution = pln.propStf.collimation.convResolution; 
                fieldWidth = pln.propStf.collimation.fieldWidth;
            else
                intConvResolution = .5; % [mm]
                fieldWidth = pln.propStf.bixelWidth;
            end

            % calculate field size and distances
            fieldLimit = ceil(fieldWidth/(2*intConvResolution));
            [F_X,F_Z] = meshgrid(-fieldLimit*intConvResolution: ...
                              intConvResolution: ...
                              (fieldLimit-1)*intConvResolution);    

            % gaussian filter to model penumbra from (measured) machine output / see
            % diploma thesis siggel 4.1.2 -> https://github.com/e0404/matRad/wiki/Dose-influence-matrix-calculation
            if isfield(obj.machine.data,'penumbraFWHMatIso')
                penumbraFWHM = obj.machine.data.penumbraFWHMatIso;
            else
                penumbraFWHM = 5;
                matRad_cfg.dispWarning('photon machine file does not contain measured penumbra width in machine.data.penumbraFWHMatIso. Assuming 5 mm.');
            end

            sigmaGauss = penumbraFWHM / sqrt(8*log(2)); % [mm] 
            % use 5 times sigma as the limits for the gaussian convolution
            gaussLimit = ceil(5*sigmaGauss/intConvResolution);
            [gaussFilterX,gaussFilterZ] = meshgrid(-gaussLimit*intConvResolution: ...
                                                intConvResolution: ...
                                                (gaussLimit-1)*intConvResolution);   
            gaussFilter =  1/(2*pi*sigmaGauss^2/intConvResolution^2) * exp(-(gaussFilterX.^2+gaussFilterZ.^2)/(2*sigmaGauss^2) );
            gaussConvSize = 2*(fieldLimit + gaussLimit);

            if ~obj.isFieldBasedDoseCalc   
                % Create fluence matrix
                F = ones(floor(fieldWidth/intConvResolution));

                if ~useCustomPrimFluenceBool
                % gaussian convolution of field to model penumbra
                    F = real(ifft2(fft2(F,gaussConvSize,gaussConvSize).*fft2(gaussFilter,gaussConvSize,gaussConvSize)));     
                end
            end

            % get kernel size and distances
            kernelLimit = ceil(lateralCutoff/intConvResolution);
            [kernelX, kernelZ] = meshgrid(-kernelLimit*intConvResolution: ...
                                        intConvResolution: ...
                                        (kernelLimit-1)*intConvResolution);

            % precalculate convoluted kernel size and distances
            kernelConvLimit = fieldLimit + gaussLimit + kernelLimit;
            [convMx_X, convMx_Z] = meshgrid(-kernelConvLimit*intConvResolution: ...
                                            intConvResolution: ...
                                            (kernelConvLimit-1)*intConvResolution);
            % calculate also the total size and distance as we need this during convolution extensively
            kernelConvSize = 2*kernelConvLimit;

            % define an effective lateral cutoff where dose will be calculated. note
            % that storage within the influence matrix may be subject to sampling
            obj.effectiveLateralCutoff = lateralCutoff + fieldWidth/2;

            counter = 0;
            matRad_cfg.dispInfo('matRad: Photon dose calculation...\n');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = 1:dij.numOfBeams % loop over all beams
                
                dij = obj.calcDoseInitBeam(ct,stf,dij,i);
                % get index of central ray or closest to the central ray
                [~,center] = min(sum(reshape([stf(i).ray.rayPos_bev],3,[]).^2));

                % get correct kernel for given SSD at central ray (nearest neighbor approximation)
                [~,currSSDIx] = min(abs([obj.machine.data.kernel.SSD]-stf(i).ray(center).SSD));
                % Display console message.
                matRad_cfg.dispInfo('\tSSD = %g mm ...\n',obj.machine.data.kernel(currSSDIx).SSD);

                kernelPos = obj.machine.data.kernelPos;
                kernel1 = obj.machine.data.kernel(currSSDIx).kernel1;
                kernel2 = obj.machine.data.kernel(currSSDIx).kernel2;
                kernel3 = obj.machine.data.kernel(currSSDIx).kernel3;

                % Evaluate kernels for all distances, interpolate between values
                kernel1Mx = interp1(kernelPos,kernel1,sqrt(kernelX.^2+kernelZ.^2),'linear',0);
                kernel2Mx = interp1(kernelPos,kernel2,sqrt(kernelX.^2+kernelZ.^2),'linear',0);
                kernel3Mx = interp1(kernelPos,kernel3,sqrt(kernelX.^2+kernelZ.^2),'linear',0);

                % convolution here if no custom primary fluence and no field based dose calc
                if ~useCustomPrimFluenceBool && ~obj.isFieldBasedDoseCalc

                    % Display console message.
                    matRad_cfg.dispInfo('\tUniform primary photon fluence -> pre-compute kernel convolution...\n');   

                    % 2D convolution of Fluence and Kernels in fourier domain
                    convMx1 = real(ifft2(fft2(F,kernelConvSize,kernelConvSize).* fft2(kernel1Mx,kernelConvSize,kernelConvSize)));
                    convMx2 = real(ifft2(fft2(F,kernelConvSize,kernelConvSize).* fft2(kernel2Mx,kernelConvSize,kernelConvSize)));
                    convMx3 = real(ifft2(fft2(F,kernelConvSize,kernelConvSize).* fft2(kernel3Mx,kernelConvSize,kernelConvSize)));

                    % Creates an interpolant for kernes from vectors position X and Z
                    if strcmp(env,'MATLAB')
                        Interp_kernel1 = griddedInterpolant(convMx_X',convMx_Z',convMx1','linear','none');
                        Interp_kernel2 = griddedInterpolant(convMx_X',convMx_Z',convMx2','linear','none');
                        Interp_kernel3 = griddedInterpolant(convMx_X',convMx_Z',convMx3','linear','none');
                    else
                        Interp_kernel1 = @(x,y)interp2(convMx_X(1,:),convMx_Z(:,1),convMx1,x,y,'linear',NaN);
                        Interp_kernel2 = @(x,y)interp2(convMx_X(1,:),convMx_Z(:,1),convMx2,x,y,'linear',NaN);
                        Interp_kernel3 = @(x,y)interp2(convMx_X(1,:),convMx_Z(:,1),convMx3,x,y,'linear',NaN);
                    end
                end

                for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!

                    counter = counter + 1;
                    obj.bixelsPerBeam = obj.bixelsPerBeam + 1;

                    % convolution here if custom primary fluence OR field based dose calc
                    if useCustomPrimFluenceBool || obj.isFieldBasedDoseCalc

                        % overwrite field opening if necessary
                        if obj.isFieldBasedDoseCalc
                            F = stf(i).ray(j).shape;
                        end

                        % prepare primary fluence array
                        primaryFluence = obj.machine.data.primaryFluence;
                        r     = sqrt( (F_X-stf(i).ray(j).rayPos(1)).^2 + (F_Z-stf(i).ray(j).rayPos(3)).^2 );
                        Psi   = interp1(primaryFluence(:,1)',primaryFluence(:,2)',r,'linear',0);

                        % apply the primary fluence to the field
                        Fx = F .* Psi;

                        % convolute with the gaussian
                        Fx = real( ifft2(fft2(Fx,gaussConvSize,gaussConvSize).* fft2(gaussFilter,gaussConvSize,gaussConvSize)) );

                        % 2D convolution of Fluence and Kernels in fourier domain
                        convMx1 = real( ifft2(fft2(Fx,kernelConvSize,kernelConvSize).* fft2(kernel1Mx,kernelConvSize,kernelConvSize)) );
                        convMx2 = real( ifft2(fft2(Fx,kernelConvSize,kernelConvSize).* fft2(kernel2Mx,kernelConvSize,kernelConvSize)) );
                        convMx3 = real( ifft2(fft2(Fx,kernelConvSize,kernelConvSize).* fft2(kernel3Mx,kernelConvSize,kernelConvSize)) );

                        % Creates an interpolant for kernes from vectors position X and Z
                        if strcmp(env,'MATLAB')
                            Interp_kernel1 = griddedInterpolant(convMx_X',convMx_Z',convMx1','linear','none');
                            Interp_kernel2 = griddedInterpolant(convMx_X',convMx_Z',convMx2','linear','none');
                            Interp_kernel3 = griddedInterpolant(convMx_X',convMx_Z',convMx3','linear','none');
                        else
                            Interp_kernel1 = @(x,y)interp2(convMx_X(1,:),convMx_Z(:,1),convMx1,x,y,'linear',NaN);
                            Interp_kernel2 = @(x,y)interp2(convMx_X(1,:),convMx_Z(:,1),convMx2,x,y,'linear',NaN);
                            Interp_kernel3 = @(x,y)interp2(convMx_X(1,:),convMx_Z(:,1),convMx3,x,y,'linear',NaN);
                        end

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
                    [ix,rad_distancesSq,isoLatDistsX,isoLatDistsZ] = matRad_calcGeoDists(obj.rot_coordsVdoseGrid, ...
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
                    bixelDose = obj.calcPhotonDoseBixel(Interp_kernel1,...
                                                               Interp_kernel2,...
                                                               Interp_kernel3,...
                                                               obj.radDepthVdoseGrid{1}(ix),...
                                                               obj.geoDistVdoseGrid{1}(ix),...
                                                               isoLatDistsX,...
                                                               isoLatDistsZ);

                    % sample dose only for bixel based dose calculation
                    if ~obj.isFieldBasedDoseCalc
                        r0   = 20 + stf(i).bixelWidth;   % [mm] sample beyond the inner core
                        Type = 'radius';
                        [ix,bixelDose] = obj.DijSampling(ix,bixelDose,obj.radDepthVdoseGrid{1}(ix),rad_distancesSq,Type,r0);
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
        
        function dij = fillDij(obj,dij,stf,pln,counter)
        % Sequentially fill the sparse matrix dij from the tmpContainer cell array
        %
        %
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
        
        function dose = calcPhotonDoseBixel(obj,Interp_kernel1,...
                  Interp_kernel2,Interp_kernel3,radDepths,geoDists,...
                  isoLatDistsX,isoLatDistsZ)
            % matRad photon dose calculation for an individual bixel
            % 
            % call
            %   dose = onj.calcPhotonDoseBixel(SAD,m,betas,Interp_kernel1,...
            %                  Interp_kernel2,Interp_kernel3,radDepths,geoDists,...
            %                  isoLatDistsX,isoLatDistsZ)
            %
            % input
            %   SAD:                source to axis distance
            %   m:                  absorption in water (part of the dose calc base
            %                       data)
            %   betas:              beta parameters for the parameterization of the 
            %                       three depth dose components
            %   Interp_kernel1/2/3: kernels for dose calculation
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

            m = obj.machine.data.m;
            betas = obj.machine.data.betas;
            % Define function_Di
            func_Di = @(beta,x) beta/(beta-m) * (exp(-m*x) - exp(-beta*x)); 

            % Calulate lateral distances using grid interpolation.
            lat1 = Interp_kernel1(isoLatDistsX,isoLatDistsZ);
            lat2 = Interp_kernel2(isoLatDistsX,isoLatDistsZ);
            lat3 = Interp_kernel3(isoLatDistsX,isoLatDistsZ);

            % now add everything together (eq 19 w/o inv sq corr -> see below)
            dose = lat1 .* func_Di(betas(1),radDepths) + ...
                   lat2 .* func_Di(betas(2),radDepths) + ...
                   lat3 .* func_Di(betas(3),radDepths);

            % inverse square correction
            dose = dose .* ((obj.machine.meta.SAD)./geoDists(:)).^2;

            % check if we have valid dose values and adjust numerical instabilities
            % from fft convolution
            dose(dose < 0 & dose > -1e-14) = 0;
            if any(isnan(dose)) || any(dose<0)
               error('Error in photon dose calculation.');
            end
        end
        
        function [ixNew,bixelDoseNew] =  DijSampling(~,ix,bixelDose,radDepthV,rad_distancesSq,sType,Param)
            % matRad dij sampling function 
            % This function samples. 
            % 
            % call
            %   [ixNew,bixelDoseNew] =
            %   matRad_DijSampling(ix,bixelDose,radDepthV,rad_distancesSq,sType,Param)
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
    
    methods (Static)

        function ret = isAvailable(pln)
            ret = any(strcmp(DoseCalcEngines.matRad_PhotonSVDPencilBeamDoseEngine.possibleRadiationModes, pln.radiationMode));
        end

    end

end

