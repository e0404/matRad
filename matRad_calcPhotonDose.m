function dij = matRad_calcPhotonDose(ct,stf,pln,cst,calcDoseDirect)
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

% initialize
matRad_calcDoseInit;

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
lateralCutoff = 50; % [mm]

% toggle custom primary fluence on/off. if 0 we assume a homogeneous
% primary fluence, if 1 we use measured radially symmetric data
useCustomPrimFluenceBool = 0;

% 0 if field calc is bixel based, 1 if dose calc is field based
isFieldBasedDoseCalc = strcmp(num2str(pln.propStf.bixelWidth),'field');

%% kernel convolution
% set up convolution grid
if isFieldBasedDoseCalc
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

% gaussian filter to model penumbra
sigmaGauss = 2.123; % [mm] / see diploma thesis siggel 4.1.2
% use 5 times sigma as the limits for the gaussian convolution
gaussLimit = ceil(5*sigmaGauss/intConvResolution);
[gaussFilterX,gaussFilterZ] = meshgrid(-gaussLimit*intConvResolution: ...
                                    intConvResolution: ...
                                    (gaussLimit-1)*intConvResolution);   
gaussFilter =  1/(2*pi*sigmaGauss^2/intConvResolution^2) * exp(-(gaussFilterX.^2+gaussFilterZ.^2)/(2*sigmaGauss^2) );
gaussConvSize = 2*(fieldLimit + gaussLimit);

if ~isFieldBasedDoseCalc   
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
effectiveLateralCutoff = lateralCutoff + fieldWidth/2;

counter = 0;
fprintf('matRad: Photon dose calculation...\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams % loop over all beams
    
    matRad_calcDoseInitBeam;
    
    % get index of central ray or closest to the central ray
    [~,center] = min(sum(reshape([stf(i).ray.rayPos_bev],3,[]).^2));
    
    % get correct kernel for given SSD at central ray (nearest neighbor approximation)
    [~,currSSDIx] = min(abs([machine.data.kernel.SSD]-stf(i).ray(center).SSD));
    
    fprintf(['                   SSD = ' num2str(machine.data.kernel(currSSDIx).SSD) 'mm                 \n']);
    
    kernelPos = machine.data.kernelPos;
    kernel1 = machine.data.kernel(currSSDIx).kernel1;
    kernel2 = machine.data.kernel(currSSDIx).kernel2;
    kernel3 = machine.data.kernel(currSSDIx).kernel3;

    % Evaluate kernels for all distances, interpolate between values
    kernel1Mx = interp1(kernelPos,kernel1,sqrt(kernelX.^2+kernelZ.^2),'linear',0);
    kernel2Mx = interp1(kernelPos,kernel2,sqrt(kernelX.^2+kernelZ.^2),'linear',0);
    kernel3Mx = interp1(kernelPos,kernel3,sqrt(kernelX.^2+kernelZ.^2),'linear',0);
    
    % convolution here if no custom primary fluence and no field based dose calc
    if ~useCustomPrimFluenceBool && ~isFieldBasedDoseCalc
        
        % Display console message.
        fprintf(['matRad: Uniform primary photon fluence -> pre-compute kernel convolution for SSD = ' ... 
                num2str(machine.data.kernel(currSSDIx).SSD) ' mm ...\n']);    

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
        bixelsPerBeam = bixelsPerBeam + 1;
    
        % convolution here if custom primary fluence OR field based dose calc
        if useCustomPrimFluenceBool || isFieldBasedDoseCalc
            
            % overwrite field opening if necessary
            if isFieldBasedDoseCalc
                F = stf(i).ray(j).shape;
            end
            
            % prepare primary fluence array
            primaryFluence = machine.data.primaryFluence;
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
        if mod(bixelsPerBeam,max(1,round(stf(i).totalNumOfBixels/200))) == 0
            matRad_progress(bixelsPerBeam/max(1,round(stf(i).totalNumOfBixels/200)),...
                            floor(stf(i).totalNumOfBixels/max(1,round(stf(i).totalNumOfBixels/200))));
        end
        % update waitbar only 100 times
        if mod(counter,round(dij.totalNumOfBixels/100)) == 0 && ishandle(figureWait)
            waitbar(counter/dij.totalNumOfBixels);
        end
        
        % remember beam and bixel number
        if ~calcDoseDirect
           dij.beamNum(counter)  = i;
           dij.rayNum(counter)   = j;
           dij.bixelNum(counter) = 1;
        else
            k = 1;
        end
        
        % Ray tracing for beam i and bixel j
        [ix,rad_distancesSq,isoLatDistsX,isoLatDistsZ] = matRad_calcGeoDists(rot_coordsVdoseGrid, ...
                                                               stf(i).sourcePoint_bev, ...
                                                               stf(i).ray(j).targetPoint_bev, ...
                                                               machine.meta.SAD, ...
                                                               find(~isnan(radDepthVdoseGrid{1})), ...
                                                               effectiveLateralCutoff);

        % empty bixels may happen during recalculation of error
        % scenarios -> skip to next bixel
        if isempty(ix)
            continue;
        end

                
        % calculate photon dose for beam i and bixel j
        bixelDose = matRad_calcPhotonDoseBixel(machine.meta.SAD,machine.data.m,...
                                                   machine.data.betas, ...
                                                   Interp_kernel1,...
                                                   Interp_kernel2,...
                                                   Interp_kernel3,...
                                                   radDepthVdoseGrid{1}(ix),...
                                                   geoDistVdoseGrid{1}(ix),...
                                                   isoLatDistsX,...
                                                   isoLatDistsZ);
                                               
        % sample dose only for bixel based dose calculation
        if ~isFieldBasedDoseCalc
            r0   = 20 + stf(i).bixelWidth;   % [mm] sample beyond the inner core
            Type = 'radius';
            [ix,bixelDose] = matRad_DijSampling(ix,bixelDose,radDepthVdoseGrid{1}(ix),rad_distancesSq,Type,r0);
        end
        
        % Save dose for every bixel in cell array
        doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(VdoseGrid(ix),1,bixelDose,dij.doseGrid.numOfVoxels,1);

        matRad_calcDoseFillDij;
               
    end
end

try
  % wait 0.1s for closing all waitbars
  allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar'); 
  delete(allWaitBarFigures);
  pause(0.1);
catch
end

