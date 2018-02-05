function dij = matRad_calcPhotonDose_4D(ct,stf,pln,cst,calcDoseDirect)
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


%{
FOR SSD FUNCTIOn
            [alpha,l,rho,~,~] = matRad_siddonRayTracer(stf(i).isoCenter, ...
                ct.resolution, ...
                stf(i).sourcePoint, ...
                stf(i).ray(j).targetPoint, ...
                ct.cube);


for k = 1:ct.numOfCtScen
    ixSSD = find(rho{k} > DensityThresholdSSD,1,'first');
    
    if isempty(ixSSD)== 1
        warning('Surface for SSD calculation starts directly in first voxel of CT\n');
    end
    
    % calculate SSD
    stf(i).ray(j).SSD{k} = 2 * stf(i).SAD * alpha(ixSSD);
end
%}

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
dij.numOfScenarios     = ct.numOfCtScen;
dij.weightToMU         = 100;
dij.scaleFactor        = 1;
dij.memorySaver        = pln.memorySaver;

% set up arrays for book keeping
dij.bixelNum = NaN*ones(dij.totalNumOfRays,1);
dij.rayNum   = NaN*ones(dij.totalNumOfRays,1);
dij.beamNum  = NaN*ones(dij.totalNumOfRays,1);

dij.nCore = cell(dij.numOfScenarios,1);
dij.nTail = cell(dij.numOfScenarios,1);
dij.nDepth = cell(dij.numOfScenarios,1);

dij.ixTail = cell(dij.numOfScenarios,1);
dij.nTailPerDepth = cell(dij.numOfScenarios,1);
dij.bixelDoseTail = cell(dij.numOfScenarios,1);

for i = 1:dij.numOfScenarios
    dij.nCore{i}    = zeros*ones(dij.totalNumOfRays,1,'uint16');
    dij.nTail{i}    = zeros*ones(dij.totalNumOfRays,1,'uint16');
    dij.nDepth{i}   = zeros*ones(dij.totalNumOfRays,1,'uint16');
    
    dij.ixTail{i}          = intmax('uint32')*ones(1000*dij.totalNumOfRays,1,'uint32');
    dij.nTailPerDepth{i}   = intmax('uint16')*ones(100*dij.totalNumOfRays,1,'uint16');
    dij.bixelDoseTail{i}   = -1*ones(100*dij.totalNumOfRays,1,'double');
end

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
V = cell(dij.numOfScenarios,1);
V{1} = [cst{:,4}];
V{1} = unique(vertcat(V{1}{:}));

% Convert CT subscripts to linear indices.
xCoordsV_vox = cell(dij.numOfScenarios,1);
yCoordsV_vox = cell(dij.numOfScenarios,1);
zCoordsV_vox = cell(dij.numOfScenarios,1);

[yCoordsV_vox{1}, xCoordsV_vox{1}, zCoordsV_vox{1}] = ind2sub(ct.cubeDim,V{1});
for i = 2:dij.numOfScenarios
    % these are probably fractional voxels
    xCoordsV_vox{i} = ct.motionVecX{i}(V{1});
    yCoordsV_vox{i} = ct.motionVecY{i}(V{1});
    zCoordsV_vox{i} = ct.motionVecZ{i}(V{1});
    
    % to get the voxel index to which each point belongs, round to the
    % nearest integer
    temp_xCoordsV_vox = round(xCoordsV_vox{i});
    temp_yCoordsV_vox = round(yCoordsV_vox{i});
    temp_zCoordsV_vox = round(zCoordsV_vox{i});
    % round up or down at boundaries
    temp_xCoordsV_vox(temp_xCoordsV_vox < 1) = 1;
    temp_xCoordsV_vox(temp_xCoordsV_vox > ct.cubeDim(2)) = ct.cubeDim(2);
    temp_yCoordsV_vox(temp_yCoordsV_vox < 1) = 1;
    temp_yCoordsV_vox(temp_yCoordsV_vox > ct.cubeDim(1)) = ct.cubeDim(1);
    temp_zCoordsV_vox(temp_zCoordsV_vox < 1) = 1;
    temp_zCoordsV_vox(temp_zCoordsV_vox > ct.cubeDim(3)) = ct.cubeDim(3);
    
    V{i} = sub2ind(ct.cubeDim,temp_yCoordsV_vox,temp_xCoordsV_vox,temp_zCoordsV_vox);
end

% set lateral cutoff value
if ~(strcmp(num2str(pln.bixelWidth),'field'))
    if calcDoseDirect
        lateralCutoff = 82; % [mm]
    else
        lateralCutoff = 50; % [mm]
    end
else
    lateralCutoff = 140; % [mm] due to the large field size
end

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

% Make a 2D grid extending +/-100mm with 0.5 mm resolution
convLimits = 100; % [mm]
convResolution = .5; % [mm]
[X,Z] = meshgrid(-convLimits:convResolution:convLimits-convResolution);

% gaussian filter to model penumbra
sigmaGauss = 2.123/convResolution; % [mm] / see diploma thesis siggel 4.1.2
gaussFilter =  1/(2*pi*sigmaGauss^2) * exp(-(X.^2+Z.^2)/(2*sigmaGauss^2*convResolution^2) );

if ~(strcmp(num2str(pln.bixelWidth),'field'))
    % Create zero matrix for the Fluence
    F = zeros(size(X));
    
    % set bixel opening to one
    F(X >= -pln.bixelWidth/2 & X < pln.bixelWidth/2 & ...
        Z >= -pln.bixelWidth/2 & Z < pln.bixelWidth/2) = 1;
    % gaussian convolution of field to model penumbra
    if ~useCustomPrimFluenceBool
        F = real(fftshift(ifft2(fft2( ifftshift(F) ).*fft2( ifftshift(gaussFilter) ))));
    end
end

counter = 0;
offsetTail = cell(dij.numOfScenarios,1);
offsetDepth = cell(dij.numOfScenarios,1);
for i = 1:dij.numOfScenarios
    offsetTail{i} = 0;
    offsetDepth{i} = 0;
end

fprintf('matRad: Photon dose calculation...\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams % loop over all beams
    
    fprintf(['Beam ' num2str(i) ' of ' num2str(dij.numOfBeams) ': \n']);
    
    bixelsPerBeam = 0;
    
    % convert voxel indices to real coordinates using iso center of beam i
    xCoordsV = cell(dij.numOfScenarios,1);
    yCoordsV = cell(dij.numOfScenarios,1);
    zCoordsV = cell(dij.numOfScenarios,1);
    coordsV = cell(dij.numOfScenarios,1);
    rot_coordsV = cell(dij.numOfScenarios,1);
    
    for j = 1:dij.numOfScenarios
        xCoordsV{j} = xCoordsV_vox{j}(:)*ct.resolution.x-stf(i).isoCenter(1);
        yCoordsV{j} = yCoordsV_vox{j}(:)*ct.resolution.y-stf(i).isoCenter(2);
        zCoordsV{j} = zCoordsV_vox{j}(:)*ct.resolution.z-stf(i).isoCenter(3);
        coordsV{j}  = [xCoordsV{j} yCoordsV{j} zCoordsV{j}];
        
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
        rot_coordsV{j} = coordsV{j}*inv_rotMx_XZ_T*inv_rotMx_XY_T;
        
        rot_coordsV{j}(:,1) = rot_coordsV{j}(:,1)-stf(i).sourcePoint_bev(1);
        rot_coordsV{j}(:,2) = rot_coordsV{j}(:,2)-stf(i).sourcePoint_bev(2);
        rot_coordsV{j}(:,3) = rot_coordsV{j}(:,3)-stf(i).sourcePoint_bev(3);
    end
    
    % ray tracing
    fprintf(['matRad: calculate radiological depth cube...']);
    [radDepthV,geoDistV] = matRad_rayTracing(stf(i),ct,V,rot_coordsV,lateralCutoff);
    fprintf('done \n');
    
    radDepthIx = cell(dij.numOfScenarios,1);
    kernel1Mx = cell(dij.numOfScenarios,1);
    kernel2Mx = cell(dij.numOfScenarios,1);
    kernel3Mx = cell(dij.numOfScenarios,1);
    Interp_kernel1 = cell(dij.numOfScenarios,1);
    Interp_kernel2 = cell(dij.numOfScenarios,1);
    Interp_kernel3 = cell(dij.numOfScenarios,1);
    for j = 1:dij.numOfScenarios
        % get indices of voxels where ray tracing results are available
        radDepthIx{j} = find(~isnan(radDepthV{j}));
        
        % limit rotated coordinates to positions where ray tracing is availabe
        rot_coordsV{j} = rot_coordsV{j}(radDepthIx{j},:);
        
        
        % get index of central ray or closest to the central ray
        [~,center] = min(sum(reshape([stf(i).ray.rayPos_bev],3,[]).^2));
        
        % get correct kernel for given SSD at central ray (nearest neighbor approximation)
        [~,currSSDIx] = min(abs([machine.data.kernel.SSD]-stf(i).ray(center).SSD{j}));
        
        fprintf(['                   SSD = ' num2str(machine.data.kernel(currSSDIx).SSD) 'mm                 \n']);
        
        kernelPos = machine.data.kernelPos;
        kernel1 = machine.data.kernel(currSSDIx).kernel1;
        kernel2 = machine.data.kernel(currSSDIx).kernel2;
        kernel3 = machine.data.kernel(currSSDIx).kernel3;
        
        % Evaluate kernels for all distances, interpolate between values
        kernel1Mx{j} = interp1(kernelPos,kernel1,sqrt(X.^2+Z.^2));
        kernel2Mx{j} = interp1(kernelPos,kernel2,sqrt(X.^2+Z.^2));
        kernel3Mx{j} = interp1(kernelPos,kernel3,sqrt(X.^2+Z.^2));
        
        % convolution here if no custom primary fluence and no field based dose calc
        if ~useCustomPrimFluenceBool && ~(strcmp(num2str(pln.bixelWidth),'field'))
            
            % Display console message.
            fprintf(['matRad: Uniform primary photon fluence -> pre-compute kernel convolution for SSD = ' ...
                num2str(machine.data.kernel(currSSDIx).SSD) ' mm ...\n']);
            
            % 2D convolution of Fluence and Kernels in fourier domain
            convMx1 = real(fftshift(ifft2(fft2( ifftshift(F) ).*fft2( ifftshift(kernel1Mx{j}) ) )));
            convMx2 = real(fftshift(ifft2(fft2( ifftshift(F) ).*fft2( ifftshift(kernel2Mx{j}) ) )));
            convMx3 = real(fftshift(ifft2(fft2( ifftshift(F) ).*fft2( ifftshift(kernel3Mx{j}) ) )));
            
            % Creates an interpolant for kernes from vectors position X and Z
            if exist('griddedInterpolant','class') % use griddedInterpoland class when available
                Interp_kernel1{j} = griddedInterpolant(X',Z',convMx1','linear');
                Interp_kernel2{j} = griddedInterpolant(X',Z',convMx2','linear');
                Interp_kernel3{j} = griddedInterpolant(X',Z',convMx3','linear');
            else
                Interp_kernel1{j} = @(x,y)interp2(X(1,:),Z(:,1),convMx1,x,y,'linear');
                Interp_kernel2{j} = @(x,y)interp2(X(1,:),Z(:,1),convMx2,x,y,'linear');
                Interp_kernel3{j} = @(x,y)interp2(X(1,:),Z(:,1),convMx3,x,y,'linear');
            end
        end
    end
    
    for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
        
        counter = counter + 1;
        bixelsPerBeam = bixelsPerBeam + 1;
        
        % convolution here if custom primary fluence OR field based dose calc
        if useCustomPrimFluenceBool || strcmp(num2str(pln.bixelWidth),'field')
            
            % get the field if field based dose calculation
            if strcmp(num2str(pln.bixelWidth),'field')
                Fx = stf(i).ray(j).shape;
                
                % apply the primary fluence to the field
            elseif useCustomPrimFluenceBool
                
                primaryFluence = machine.data.primaryFluence;
                r     = sqrt( (X-stf(i).ray(j).rayPos(1)).^2 + (Z-stf(i).ray(j).rayPos(3)).^2 );
                Psi   = interp1(primaryFluence(:,1)',primaryFluence(:,2)',r);
                Fx = F .* Psi;
                
            end
            
            % convolute with the gaussian
            Fx = real(fftshift(ifft2(fft2( ifftshift(Fx) ).*fft2( ifftshift(gaussFilter) ))));
            
            for k = 1:dij.numOfScenarios
                % 2D convolution of Fluence and Kernels in fourier domain
                convMx1 = real(fftshift(ifft2(fft2( ifftshift(Fx) ).*fft2( ifftshift(kernel1Mx{k}) ) )));
                convMx2 = real(fftshift(ifft2(fft2( ifftshift(Fx) ).*fft2( ifftshift(kernel2Mx{k}) ) )));
                convMx3 = real(fftshift(ifft2(fft2( ifftshift(Fx) ).*fft2( ifftshift(kernel3Mx{k}) ) )));
                
                % Creates an interpolant for kernes from vectors position X and Z
                if exist('griddedInterpolant','class') % use griddedInterpoland class when available
                    Interp_kernel1{k} = griddedInterpolant(X',Z',convMx1','linear');
                    Interp_kernel2{k} = griddedInterpolant(X',Z',convMx2','linear');
                    Interp_kernel3{k} = griddedInterpolant(X',Z',convMx3','linear');
                else
                    Interp_kernel1{k} = @(x,y)interp2(X(1,:),Z(:,1),convMx1,x,y,'linear');
                    Interp_kernel2{k} = @(x,y)interp2(X(1,:),Z(:,1),convMx2,x,y,'linear');
                    Interp_kernel3{k} = @(x,y)interp2(X(1,:),Z(:,1),convMx3,x,y,'linear');
                end
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
        
        for k = 1:dij.numOfScenarios
            % calculate photon dose for beam i and bixel j
            bixelDose = matRad_calcPhotonDoseBixel(machine.meta.SAD,machine.data.m,...
                machine.data.betas, ...
                Interp_kernel1{k},...
                Interp_kernel2{k},...
                Interp_kernel3{k},...
                radDepthV{k}(ix{k}),...
                geoDistV{k}(ix{k}),...
                isoLatDistsX{k},...
                isoLatDistsZ{k});
            
            % sample dose only for bixel based dose calculation, and only if
            % not calculating dose directly!
            if ~strcmp(num2str(pln.bixelWidth),'field') && ~calcDoseDirect
                r0   = 25;   % [mm] sample beyond the iresnner core
                Type = 'radius';
                
                if dij.memorySaver
                    [ix{k},bixelDose,ixTail,nTailPerDepth,bixelDoseTail,nTail,nDepth,nCore] = matRad_DijSampling_memorySaver(ix{k},bixelDose,radDepthV{k}(ix{k}),rad_distancesSq{k},Type,r0);
                    
                    dij.ixTail{k}(offsetTail{k}+(1:nTail)) = uint32(V{1}(ixTail));
                    dij.nTailPerDepth{k}(offsetDepth{k}+(1:nDepth)) = nTailPerDepth;
                    dij.bixelDoseTail{k}(offsetDepth{k}+(1:nDepth)) = bixelDoseTail;
                    
                    dij.nCore{k}(counter) = uint16(nCore);
                    dij.nTail{k}(counter) = uint16(nTail);
                    dij.nDepth{k}(counter) = uint16(nDepth);
                    
                    offsetTail{k} = offsetTail{k}+nTail;
                    offsetDepth{k} = offsetDepth{k}+nDepth;
                else
                    [ix{k},bixelDose] = matRad_DijSampling(ix{k},bixelDose,radDepthV{k}(ix{k}),rad_distancesSq{k},Type,r0);
                end
            end
            % Save dose for every bixel in cell array
            doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,k} = sparse(V{1}(ix{k}),1,bixelDose,dij.numOfVoxels,1);
            
            % save computation time and memory by sequentially filling the
            % sparse matrix dose.dij from the cell array
            if mod(counter,numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
                if calcDoseDirect
                    if isfield(stf(1).ray(1),'weight')
                        % score physical dose
                        dij.physicalDose{k}(:,1) = dij.physicalDose{k}(:,1) + stf(i).ray(j).weight * doseTmpContainer{1,k};
                    else
                        error(['No weight available for beam ' num2str(i) ', ray ' num2str(j)]);
                    end
                else
                    % fill entire dose influence matrix
                    dij.physicalDose{k}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [doseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,k}];
                end
            end
        end
        
    end
end

for i = 1:dij.numOfScenarios
    dij.ixTail{i}(dij.ixTail{i} == intmax('uint32')) = [];
    dij.nTailPerDepth{i}(dij.nTailPerDepth{i} == intmax('uint16')) = [];
    dij.bixelDoseTail{i}(dij.bixelDoseTail{i} == -1) = [];
end

try
    % wait 0.1s for closing all waitbars
    allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar');
    delete(allWaitBarFigures);
    pause(0.1);
catch
end


