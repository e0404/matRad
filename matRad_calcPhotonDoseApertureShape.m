function resultGUI = matRad_calcPhotonDoseApertureShape(ct,cst,pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad aperture based dose calcualtion
% 
% call
%   d = matRad_calcPhotonDoseApertureShape(ct,cst,pln,visMode)
%
% input
%
%
% output
%
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    visMode = 0;
end

lateralCutOff = 100; % [mm] / calculate dose within a corridor arround the target
rayMxSpacing  = min(ct.resolution);
cornerPoints  = [                     1                      1                      1;
                                      1                      1 pln.voxelDimensions(3);
                                      1 pln.voxelDimensions(2) pln.voxelDimensions(3);
                                      1 pln.voxelDimensions(2)                      1;
                 pln.voxelDimensions(1)                      1                      1;
                 pln.voxelDimensions(1)                      1 pln.voxelDimensions(3);
                 pln.voxelDimensions(1) pln.voxelDimensions(2) pln.voxelDimensions(3);
                 pln.voxelDimensions(1) pln.voxelDimensions(2)                      1];

distToCornerPoints = sqrt(sum(((cornerPoints-repmat(pln.isoCenter./ct.resolution,8,1)+1).*repmat(ct.resolution,8,1)).^2,2));

rayMxDistBehindIsocenterPlane = max(distToCornerPoints);
             
% throw error for particles
if strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    error('aperture based dose calculation only possible for photons\n');
elseif strcmp(pln.radiationMode,'photons')
    load photonPencilBeamKernels_6MV;
end

% find all target voxels from cst cell array
targetVoxels = [];
for i=1:size(cst,1)
    if isequal(cst{i,3},'TARGET')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% attention hack think if this should really be neglected at later point && ~isempty(cst{i,6})
        targetVoxels = [targetVoxels;cst{i,4}];
    end
end

% Remove double voxels
targetVoxels = unique(targetVoxels);

% generate binary target cube
voi    = zeros(size(ct.cube));
voi(targetVoxels) = 1;
    
% add margin
addmarginBool = 1;
if addmarginBool
    voi          = matRad_addMargin(voi,ct.resolution,lateralCutOff*[1 1 1],true);
    targetVoxels = find(voi>0);
end

% throw error message if no target is found
if isempty(targetVoxels)
    error('Could not find target');
end

% Convert linear indices to 3D voxel coordinates
[coordsY_targetVoxels, coordsX_targetVoxels, coordsZ_targetVoxels] = ind2sub(size(ct.cube),targetVoxels);

coordsX_targetVoxels = coordsX_targetVoxels(:)*ct.resolution(1) - pln.isoCenter(1);
coordsY_targetVoxels = coordsY_targetVoxels(:)*ct.resolution(2) - pln.isoCenter(2);
coordsZ_targetVoxels = coordsZ_targetVoxels(:)*ct.resolution(3) - pln.isoCenter(3);
coords_targetVoxels  = [coordsX_targetVoxels coordsY_targetVoxels coordsZ_targetVoxels];

% voxels inside patient
patientVoxels = unique([cell2mat(cst(:,4))]);

% Convert CT subscripts to linear indices.
[coordsY_patientVoxels, coordsX_patientVoxels, coordsZ_patientVoxels] = ind2sub(size(ct.cube),patientVoxels);

coordsX_patientVoxels = coordsX_patientVoxels*ct.resolution(1) - pln.isoCenter(1);
coordsY_patientVoxels = coordsY_patientVoxels*ct.resolution(2) - pln.isoCenter(2);
coordsZ_patientVoxels = coordsZ_patientVoxels*ct.resolution(3) - pln.isoCenter(3);
coords_patientVoxels  = [coordsX_patientVoxels coordsY_patientVoxels coordsZ_patientVoxels];

% all voxels
[X_geo,Y_geo,Z_geo] = meshgrid(ct.resolution(1)*(1:size(ct.cube,2)),...
    ct.resolution(2)*(1:size(ct.cube,1)),ct.resolution(3)*(1:size(ct.cube,3)));

coords_all = [X_geo(:)-pln.isoCenter(1) Y_geo(:)-pln.isoCenter(2) Z_geo(:)-pln.isoCenter(3)];

% other results
geoDistCube      = NaN*ones(pln.voxelDimensions);
latDistXCube     = NaN*ones(pln.voxelDimensions);
latDistZCube     = NaN*ones(pln.voxelDimensions);
radDepthCube     = NaN*ones(pln.voxelDimensions);

% Define steering file like struct. Prellocating for speed.
stf = struct;

% loop over all angles
for i = 1:length(pln.gantryAngles)
    
    % Save meta information for treatment plan
    stf(i).gantryAngle   = pln.gantryAngles(i);
    stf(i).couchAngle    = pln.couchAngles(i);
    stf(i).bixelWidth    = pln.bixelWidth;
    stf(i).radiationMode = pln.radiationMode;
    
    % gantry and couch roation matrices according to IEC 61217 standard
    % instead of moving the beam around the patient, we perform an inverse
    % rotation of the patient, i.e. we consider a beam's eye view
    % coordinate system
    
    % Rotation around Z axis (gantry)
    rotMx_XY = [cosd(pln.gantryAngles(i)) -sind(pln.gantryAngles(i)) 0;
                sind(pln.gantryAngles(i))  cosd(pln.gantryAngles(i)) 0;
                                        0                          0 1];
    
    % Rotation around Y axis (Couch movement)
    rotMx_XZ = [ cosd(pln.couchAngles(i)) 0 sind(pln.couchAngles(i));
                                        0 1                         0;
                -sind(pln.couchAngles(i)) 0  cosd(pln.couchAngles(i))];
    
    % rotate target coordinates around Y axis and then around Z axis
    % i.e. 1st couch, 2nd gantry; matrix multiplication not cummutative
    rot_coords_targetVoxels = coords_targetVoxels*rotMx_XZ*rotMx_XY;
    
    % project x and z coordinates to isocenter
    coordsAtIsoCenterPlane_targetVoxels(:,1) = (rot_coords_targetVoxels(:,1)*pln.SAD)./(pln.SAD + rot_coords_targetVoxels(:,2));
    coordsAtIsoCenterPlane_targetVoxels(:,2) = (rot_coords_targetVoxels(:,3)*pln.SAD)./(pln.SAD + rot_coords_targetVoxels(:,2));
    
    % Take unique rows values for beamlets positions. Calculate position of
    % central ray for every bixel    
    rayPos = unique(rayMxSpacing*round([            coordsAtIsoCenterPlane_targetVoxels(:,1) ... 
                                         zeros(size(coordsAtIsoCenterPlane_targetVoxels,1),1) ...
                                                    coordsAtIsoCenterPlane_targetVoxels(:,2)]/rayMxSpacing),'rows');
                                                
    % pad ray position array if resolution of target voxel grid not sufficient                                                
	if rayMxSpacing<max(ct.resolution)
        origRayPos = rayPos;
        for j = -floor(max(ct.resolution)/rayMxSpacing):floor(max(ct.resolution)/rayMxSpacing)
            for k = -floor(max(ct.resolution)/rayMxSpacing):floor(max(ct.resolution)/rayMxSpacing)
                if abs(j)+abs(k)==0
                    continue;
                end
                
                rayPos = [rayPos; origRayPos(:,1)+j*rayMxSpacing origRayPos(:,2) origRayPos(:,3)+k*rayMxSpacing];
                                
            end
        end
    end
    
    % remove double rays
    rayPos = unique(rayPos,'rows');
    
    % Save the number of rays
    stf(i).numOfRays = size(rayPos,1);
    
    % Save ray and target position in beam eye´s view (bev)
    for j = 1:stf(i).numOfRays
        stf(i).ray(j).rayPos_bev = rayPos(j,:);
        stf(i).ray(j).targetPoint_bev = [stf(i).ray(j).rayPos_bev(1) ...
                                         rayMxDistBehindIsocenterPlane ...
                                         stf(i).ray(j).rayPos_bev(3)];
    end
    
    % source position in bev
    sourcePoint_bev = [0 -pln.SAD 0];
    
    % compute coordinates in lps coordinate system, i.e. rotate beam
    % geometry around fixed patient
    
    % Rotation around Z axis (gantry)
    rotMx_XY_rotated = [ cosd(pln.gantryAngles(i)) sind(pln.gantryAngles(i)) 0;
                        -sind(pln.gantryAngles(i)) cosd(pln.gantryAngles(i)) 0;
                                                 0                         0 1];
    
    % Rotation around Y axis (couch)
    rotMx_XZ_rotated = [ cosd(pln.couchAngles(i)) 0 -sind(pln.couchAngles(i));
                                                0 1                        0;
                         sind(pln.couchAngles(i)) 0 cosd(pln.couchAngles(i))];
    
    % Rotated Source point, first needs to be rotated around gantry, and then
    % couch.
    stf(i).sourcePoint = sourcePoint_bev*rotMx_XY_rotated*rotMx_XZ_rotated;
    
    % Save ray and target position in lps system.
    for j = 1:stf(i).numOfRays
        stf(i).ray(j).rayPos      = stf(i).ray(j).rayPos_bev*rotMx_XY_rotated*rotMx_XZ_rotated;
        stf(i).ray(j).targetPoint = stf(i).ray(j).targetPoint_bev*rotMx_XY_rotated*rotMx_XZ_rotated;
    end
    
    % now rotate all patient voxels in bev to determine corresponding rays
    rot_coords_patientVoxels = coords_patientVoxels*rotMx_XZ*rotMx_XY;
    
    rot_coords_patientVoxels(:,1) = rot_coords_patientVoxels(:,1)-sourcePoint_bev(1);
    rot_coords_patientVoxels(:,2) = rot_coords_patientVoxels(:,2)-sourcePoint_bev(2);
    rot_coords_patientVoxels(:,3) = rot_coords_patientVoxels(:,3)-sourcePoint_bev(3);
    
    
    % calculate lateral distances and geometric distances for central ray
    [ix,~,geoDists,x_latDists,z_latDists] = matRad_calcRadGeoDists(ct.cube, ...
                                                patientVoxels, ...
                                                pln.isoCenter, ...
                                                rot_coords_patientVoxels, ...
                                                ct.resolution, ...
                                                stf(i).sourcePoint, ...
                                                [0 rayMxDistBehindIsocenterPlane 0]*rotMx_XY_rotated*rotMx_XZ_rotated, ...
                                                sourcePoint_bev, ...
                                                [0 rayMxDistBehindIsocenterPlane 0], ...
                                                coords_patientVoxels, ...
                                                inf);
        
    latDistXCube(patientVoxels(ix)) = x_latDists;
    latDistZCube(patientVoxels(ix)) = z_latDists;
    geoDistCube (patientVoxels(ix)) = geoDists;

    % calc radiological depths with ray tracing of entire cube
    for j = stf(i).numOfRays:-1:1
    
        stf(i).ray(j).energy = energy;
        stf(i).numOfBixelsPerRay(j) = NaN;
        
        [ix,radDepths] = matRad_calcRadDistsRay(ct.cube, ...
                            pln.isoCenter, ...
                            ct.resolution, ...
                            stf(i).sourcePoint, ...
                            stf(i).ray(j).targetPoint, ...
                            coords_all);
                                    
        radDepthCube(ix) = radDepths;
        
    end
        
    %{
    for i = 1:size(ct.cube,3)
        clf
        hold on
        subplot(2,2,1)
        imagesc(squeeze(radDepthCube(:,:,i)));
        colorbar
        axis equal tight
        title(['radDepth ' num2str(i)])
        subplot(2,2,2)
        imagesc(squeeze(geoDistCube(:,:,i)));
        colorbar
        axis equal tight
        title(['geoDist ' num2str(i)])
        subplot(2,2,3)
        imagesc(squeeze(latDistZCube(:,:,i)));
        colorbar
        axis equal tight
        title(['latDistZ ' num2str(i)])
        subplot(2,2,4)
        imagesc(squeeze(latDistXCube(:,:,i)));
        colorbar
        axis equal tight
        title(['latDistX ' num2str(i)])
        drawnow
        pause
    end
    %}
    
    % save total number of bixels
    stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % actual dose calculation
       
    % Make a 2D grid extending +/-100mm with 0.1 mm resolution
    convLimits = 100; % [mm]
    convResolution = .5; % [mm]
    [X,Z] = meshgrid(-convLimits:convResolution:convLimits);

    % Create zero matrix for the Fluence
    F = zeros(size(X));

    % set shape
    fieldSize = 50;
    F(abs(X)<=fieldSize & abs(Z)<=fieldSize) = 1;
    
    % gaussian convolution of field to model penumbra
    sigmaGauss = 2.1/convResolution; % [mm] / see diploma thesis siggel 4.1.2
    gaussFilter =  convResolution^2/(2*pi*sigmaGauss^2) * exp( -(X.^2+Z.^2)/(2*sigmaGauss^2) );
    F = real(fftshift(ifft2(fft2( ifftshift(F) ).*fft2( ifftshift(gaussFilter) ))));

    % multiplication with primary fluence
    r     = sqrt( X.^2 + Z.^2 );
    Psi   = interp1(primaryFluence(:,1),primaryFluence(:,2),r);
    FxPsi = F .* Psi;
        
    % Evaluate piecewise polynomial kernels
    kernel1Mx = ppval(ppKernel1,sqrt(X.^2+Z.^2));
    kernel2Mx = ppval(ppKernel2,sqrt(X.^2+Z.^2));
    kernel3Mx = ppval(ppKernel3,sqrt(X.^2+Z.^2));

    % 2D convolution of Fluence and Kernels in fourier domain
    convMx1 = real(fftshift(ifft2(fft2( ifftshift(FxPsi) ).*fft2( ifftshift(kernel1Mx) ) )));
    convMx2 = real(fftshift(ifft2(fft2( ifftshift(FxPsi) ).*fft2( ifftshift(kernel2Mx) ) )));
    convMx3 = real(fftshift(ifft2(fft2( ifftshift(FxPsi) ).*fft2( ifftshift(kernel3Mx) ) )));
    
    % Creates an interpolant for kernes from vectors position X and Z
    if exist('griddedInterpolant','class') % use griddedInterpoland class when available 
        Interp_kernel1 = griddedInterpolant(X',Z',convMx1','linear');
        Interp_kernel2 = griddedInterpolant(X',Z',convMx2','linear');
        Interp_kernel3 = griddedInterpolant(X',Z',convMx3','linear');
    else
        Interp_kernel1 = @(x,y)interp2(X(1,:),Z(:,1),convMx1,x,y,'linear');
        Interp_kernel2 = @(x,y)interp2(X(1,:),Z(:,1),convMx2,x,y,'linear');
        Interp_kernel3 = @(x,y)interp2(X(1,:),Z(:,1),convMx3,x,y,'linear');
    end
    
    % scale lateral distances to iso center plane
    latDistXCube = (latDistXCube) ./ geoDistCube .* pln.SAD;
    latDistZCube = (latDistZCube) ./ geoDistCube .* pln.SAD;
       
    % Calulate lateral distances using grid interpolation.
    lat1 = Interp_kernel1(latDistXCube,latDistZCube);
    lat2 = Interp_kernel2(latDistXCube,latDistZCube);
    lat3 = Interp_kernel3(latDistXCube,latDistZCube);
    
    % Define function_Di
    func_Di = @(beta,x) beta/(beta-m) * (exp(-m*x) - exp(-beta*x)); 

    % now add everything together (eq 19 w/o inv sq corr -> see below)
    dose = lat1 .* func_Di(betas(1),radDepthCube) + ...
           lat2 .* func_Di(betas(2),radDepthCube) + ...
           lat3 .* func_Di(betas(3),radDepthCube);

    % inverse square correction
    resultGUI.physicalDose = dose .* (pln.SAD./geoDistCube).^2;
       
    % Show progress
    matRad_progress(i,length(pln.gantryAngles));
    
end