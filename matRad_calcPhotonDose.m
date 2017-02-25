function dij = matRad_calcPhotonDose(ct,stf,pln,cst,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% set consistent random seed (enables reproducibility)
rng(0);

% default: dose influence matrix computation
if ~exist('calcDoseDirect','var')
    calcDoseDirect = false;
end

% issue warning if biological optimization not possible
if sum(strcmp(pln.bioOptimization,{'effect','RBExD'}))>0
    warndlg('Effect based and RBE optimization not available for photons - physical optimization is carried out instead.');
    pln.bioOptimization = 'none';
end

% initialize waitbar
figureWait = waitbar(0,'calculate dose influence matrix for photons...');
% show busy state
set(figureWait,'pointer','watch');

% meta information for dij
dij.numOfBeams         = pln.numOfBeams;
dij.numOfVoxels        = pln.numOfVoxels;
dij.resolution         = ct.resolution;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.dimensions         = pln.voxelDimensions;
dij.numOfScenarios     = 1;

% set up arrays for book keeping
dij.bixelNum = NaN*ones(dij.totalNumOfRays,1);
dij.rayNum   = NaN*ones(dij.totalNumOfRays,1);
dij.beamNum  = NaN*ones(dij.totalNumOfRays,1);

% Allocate space for dij.physicalDose sparse matrix
for i = 1:dij.numOfScenarios
    dij.physicalDose{i} = spalloc(prod(ct.cubeDim),dij.totalNumOfBixels,1);
end

% Allocate memory for dose_temp cell array
if calcDoseDirect
    numOfBixelsContainer = 1;
else
    numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
end
doseTmpContainer = cell(numOfBixelsContainer,dij.numOfScenarios);

% take only voxels inside patient
V = [cst{:,4}];
V = unique(vertcat(V{:}));

% ignore densities outside of contours
eraseCtDensMask = ones(dij.numOfVoxels,1);
eraseCtDensMask(V) = 0;
for i = 1:dij.numOfScenarios
    ct.cube{i}(eraseCtDensMask == 1) = 0;
end

% Convert CT subscripts to linear indices.
[yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(ct.cubeDim,V);

% set lateral cutoff value
lateralCutoff = 82; % [mm]

% toggle custom primary fluence on/off. if 0 we assume a homogeneous
% primary fluence, if 1 we use measured radially symmetric data
useCustomPrimFluenceBool = 1;

%% kernel convolution
% prepare data for convolution to reduce calculation time
fileName = [pln.radiationMode '_' pln.machine];
try
   load([fileparts(mfilename('fullpath')) filesep fileName]);
catch
   error(['Could not find the following machine file: ' fileName ]); 
end

% set up convolution grid
if strcmp(num2str(pln.bixelWidth),'field')
    % get data from DICOM import
    intConvResolution = pln.Collimation.convResolution; 
    fieldWidth = pln.Collimation.fieldWidth;
else
    intConvResolution = .5; % [mm]
    fieldWidth = pln.bixelWidth;
end

% Make a 2D grid extending +/- lateral cutoff for kernels
intConvLimits = lateralCutoff; % [mm]
[kernelX,kernelZ] = meshgrid(-intConvLimits:intConvResolution:intConvLimits-intConvResolution);   

% gaussian filter to model penumbra
sigmaGauss = 2.123/intConvResolution; % [mm] / see diploma thesis siggel 4.1.2
gaussFilter =  1/(2*pi*sigmaGauss^2) * exp(-(kernelX.^2+kernelZ.^2)/(2*sigmaGauss^2*intConvResolution^2) );
kernelSize = size(gaussFilter,1);

if ~(strcmp(num2str(pln.bixelWidth),'field'))
    
    % Create fluence matrix
    F = ones(floor(pln.bixelWidth/intConvResolution));
    Fsize = size(F,1);
    
    if ~useCustomPrimFluenceBool
    % gaussian convolution of field to model penumbra
        F = ifft2(fft2(F,Fsize+kernelSize,Fsize+kernelSize).* fft2(gaussFilter,Fsize+kernelSize,Fsize+kernelSize));
        Fsize = size(F,1);
    
        convOffset = intConvResolution*ceil((Fsize-1)./2);  
    
        [convMx_X,convMx_Z] = meshgrid(-(intConvLimits+convOffset): ...
                         intConvResolution: ...
                         (intConvLimits-intConvResolution+convOffset));
    else
        convOffset = intConvResolution*ceil((Fsize-1)./2);  
    
        [F_X,F_Z] = meshgrid(-convOffset: ...
                         intConvResolution: ...
                         convOffset-intConvResolution);    
    end
end

counter = 0;

fprintf('matRad: Photon dose calculation...\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams % loop over all beams
    
    fprintf(['Beam ' num2str(i) ' of ' num2str(dij.numOfBeams) ': \n']);

    bixelsPerBeam = 0;

    % convert voxel indices to real coordinates using iso center of beam i
    xCoordsV = xCoordsV_vox(:)*ct.resolution.x-stf(i).isoCenter(1);
    yCoordsV = yCoordsV_vox(:)*ct.resolution.y-stf(i).isoCenter(2);
    zCoordsV = zCoordsV_vox(:)*ct.resolution.z-stf(i).isoCenter(3);
    coordsV  = [xCoordsV yCoordsV zCoordsV];

    % Set gantry and couch rotation matrices according to IEC 61217
    % Use transpose matrices because we are working with row vectros

    % rotation around Z axis (gantry)
    inv_rotMx_XY_T = [ cosd(-pln.gantryAngles(i)) sind(-pln.gantryAngles(i)) 0;
                      -sind(-pln.gantryAngles(i)) cosd(-pln.gantryAngles(i)) 0;
                                                0                          0 1];

    % rotation around Y axis (couch)
    inv_rotMx_XZ_T = [cosd(-pln.couchAngles(i)) 0 -sind(-pln.couchAngles(i));
                                              0 1                         0;
                      sind(-pln.couchAngles(i)) 0  cosd(-pln.couchAngles(i))];

    % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
    rot_coordsV = coordsV*inv_rotMx_XZ_T*inv_rotMx_XY_T;

    rot_coordsV(:,1) = rot_coordsV(:,1)-stf(i).sourcePoint_bev(1);
    rot_coordsV(:,2) = rot_coordsV(:,2)-stf(i).sourcePoint_bev(2);
    rot_coordsV(:,3) = rot_coordsV(:,3)-stf(i).sourcePoint_bev(3);

    % ray tracing
    fprintf('matRad: calculate radiological depth cube...');
    [radDepthV,geoDistV] = matRad_rayTracing(stf(i),ct,V,rot_coordsV,intConvLimits);
    fprintf('done \n');
    
    % get indices of voxels where ray tracing results are available
    radDepthIx = find(~isnan(radDepthV{1}));
    
    % limit rotated coordinates to positions where ray tracing is availabe
    rot_coordsV = rot_coordsV(radDepthIx,:);
    
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
    if ~useCustomPrimFluenceBool && ~(strcmp(num2str(pln.bixelWidth),'field'))
        
        % Display console message.
        fprintf(['matRad: Uniform primary photon fluence -> pre-compute kernel convolution for SSD = ' ... 
                num2str(machine.data.kernel(currSSDIx).SSD) ' mm ...\n']);    

        % 2D convolution of Fluence and Kernels in fourier domain
        convMx1 = ifft2(fft2(F,Fsize+kernelSize,Fsize+kernelSize).* fft2(kernel1Mx,Fsize+kernelSize,Fsize+kernelSize));
        convMx2 = ifft2(fft2(F,Fsize+kernelSize,Fsize+kernelSize).* fft2(kernel2Mx,Fsize+kernelSize,Fsize+kernelSize));
        convMx3 = ifft2(fft2(F,Fsize+kernelSize,Fsize+kernelSize).* fft2(kernel3Mx,Fsize+kernelSize,Fsize+kernelSize));

        % Creates an interpolant for kernes from vectors position X and Z
        if exist('griddedInterpolant','class') % use griddedInterpoland class when available 
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
        if useCustomPrimFluenceBool || strcmp(num2str(pln.bixelWidth),'field')
            
            % overwrite field opening if necessary
            if strcmp(num2str(pln.bixelWidth),'field')
                F = stf(i).ray(j).shape;
                Fsize = size(F,1);
                convOffset = intConvResolution*ceil((Fsize-1)./2);  
    
                [F_X,F_Z] = meshgrid(-convOffset: ...
                                 intConvResolution: ...
                                 convOffset-intConvResolution);
            end
            
            % prepare primary fluence array
            primaryFluence = machine.data.primaryFluence;
            r     = sqrt( (F_X-stf(i).ray(j).rayPos(1)).^2 + (F_Z-stf(i).ray(j).rayPos(3)).^2 );
            Psi   = interp1(primaryFluence(:,1)',primaryFluence(:,2)',r,'linear',0);
                
            % apply the primary fluence to the field
            Fx = F .* Psi;
            
            % convolute with the gaussian
            Fx = real( ifft2(fft2(Fx,Fsize+kernelSize,Fsize+kernelSize).* fft2(gaussFilter,Fsize+kernelSize,Fsize+kernelSize)) );
            Fxsize = size(Fx,1);

            % 2D convolution of Fluence and Kernels in fourier domain
            convMx1 = real( ifft2(fft2(Fx,Fxsize+kernelSize,Fxsize+kernelSize).* fft2(kernel1Mx,Fxsize+kernelSize,Fxsize+kernelSize)) );
            convMx2 = real( ifft2(fft2(Fx,Fxsize+kernelSize,Fxsize+kernelSize).* fft2(kernel2Mx,Fxsize+kernelSize,Fxsize+kernelSize)) );
            convMx3 = real( ifft2(fft2(Fx,Fxsize+kernelSize,Fxsize+kernelSize).* fft2(kernel3Mx,Fxsize+kernelSize,Fxsize+kernelSize)) );

            % set up grid
            convOffset = intConvResolution*ceil((Fxsize-1)./2);  
            [convMx_X,convMx_Z] = meshgrid(-(intConvLimits+convOffset): ...
                             intConvResolution: ...
                             (intConvLimits-intConvResolution+convOffset));

            % Creates an interpolant for kernes from vectors position X and Z
            if exist('griddedInterpolant','class') % use griddedInterpoland class when available 
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
        dij.beamNum(counter)  = i;
        dij.rayNum(counter)   = j;
        dij.bixelNum(counter) = j;
        
        % Ray tracing for beam i and bixel j
        [ix,rad_distancesSq,isoLatDistsX,isoLatDistsZ] = matRad_calcGeoDists(rot_coordsV, ...
                                                               stf(i).sourcePoint_bev, ...
                                                               stf(i).ray(j).targetPoint_bev, ...
                                                               machine.meta.SAD, ...
                                                               radDepthIx, ...
                                                               lateralCutoff);

        % calculate photon dose for beam i and bixel j
        bixelDose = matRad_calcPhotonDoseBixel(machine.meta.SAD,machine.data.m,...
                                                   machine.data.betas, ...
                                                   Interp_kernel1,...
                                                   Interp_kernel2,...
                                                   Interp_kernel3,...
                                                   radDepthV{1}(ix),...
                                                   geoDistV(ix),...
                                                   isoLatDistsX,...
                                                   isoLatDistsZ);
                                               
        % sample dose only for bixel based dose calculation
        if ~strcmp(num2str(pln.bixelWidth),'field')
            r0   = 25;   % [mm] sample beyond the inner core
            Type = 'radius';
            [ix,bixelDose] = matRad_DijSampling(ix,bixelDose,radDepthV{1}(ix),rad_distancesSq,Type,r0);
        end   
        % Save dose for every bixel in cell array
        doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(V(ix),1,bixelDose,dij.numOfVoxels,1);
                
        % save computation time and memory by sequentially filling the 
        % sparse matrix dose.dij from the cell array
        if mod(counter,numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
            if calcDoseDirect
                if isfield(stf(1).ray(1),'weight')
                    % score physical dose
                    dij.physicalDose{1}(:,1) = dij.physicalDose{1}(:,1) + stf(i).ray(j).weight * doseTmpContainer{1,1};
                else
                    error(['No weight available for beam ' num2str(i) ', ray ' num2str(j)]);
                end
            else
                % fill entire dose influence matrix
                dij.physicalDose{1}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [doseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
            end
        end
        
    end
end

try
  % wait 0.1s for closing all waitbars
  allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar'); 
  delete(allWaitBarFigures);
  pause(0.1);
catch
end


