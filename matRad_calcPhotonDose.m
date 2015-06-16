function dij = matRad_calcPhotonDose(ct,stf,pln,cst,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad photon dose calculation wrapper
% 
% call
%   dij = matRad_calcPhotonDose(ct,stf,pln,cst,visBool)
%
% input
%   ct:         ct cube
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%   cst:        matRad cst struct
%   visBool:    toggle on/off visualization (optional)
%
% output
%   dij:        matRad dij struct
%
% References
%   [1] http://www.ncbi.nlm.nih.gov/pubmed/8497215
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

% if visBool not set toogle off visualization
if nargin < 5
    visBool = 0;
end

% initialize waitbar
figureWait=waitbar(0,'photon dij-calculation..');
% meta information for dij
dij.numOfBeams         = pln.numOfBeams;
dij.numOfVoxels        = pln.numOfVoxels;
dij.resolution         = ct.resolution;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.dimensions         = pln.voxelDimensions;

% set up arrays for book keeping
dij.bixelNum = NaN*ones(dij.totalNumOfRays,1);
dij.rayNum   = NaN*ones(dij.totalNumOfRays,1);
dij.beamNum  = NaN*ones(dij.totalNumOfRays,1);

% Allocate space for dij.physicalDose sparse matrix
dij.physicalDose = spalloc(numel(ct.cube),dij.totalNumOfBixels,1);

% Allocate memory for dose_temp cell array
numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
doseTmpContainer = cell(numOfBixelsContainer,1);

% take only voxels inside patient
V = unique([cell2mat(cst(:,4))]);

% Convert CT subscripts to linear indices.
[yCoordsV, xCoordsV, zCoordsV] = ind2sub(size(ct.cube),V);

xCoordsV = (xCoordsV(:)-0.5)*ct.resolution(1)-pln.isoCenter(1);
yCoordsV = (yCoordsV(:)-0.5)*ct.resolution(2)-pln.isoCenter(2);
zCoordsV = (zCoordsV(:)-0.5)*ct.resolution(3)-pln.isoCenter(3);
coords_inside = [xCoordsV yCoordsV zCoordsV];

% set lateral cutoff value
lateralCutoff = 20; % [mm]

%% kernel convolution
% load polynomial fits for kernels ppKernel1, ppKernel2, ppKernel3
load photonPencilBeamKernels_6MV.mat;

% Display console message.
fprintf('matRad: Kernel convolution... \n');

% Make a 2D grid extending +/-100mm with 0.1 mm resolution
[X,Z] = meshgrid(-100:0.1:100);

% Evaluate piecewise polynomial kernels
kernel1Mx = ppval(ppKernel1,sqrt(X.^2+Z.^2));
kernel2Mx = ppval(ppKernel2,sqrt(X.^2+Z.^2));
kernel3Mx = ppval(ppKernel3,sqrt(X.^2+Z.^2));

% Create zero matrix for the Fluence
F = zeros(size(X));

% set bixel opening to one
F(1001-pln.bixelWidth/2/0.1:1001+pln.bixelWidth/2/0.1,...
  1001-pln.bixelWidth/2/0.1:1001+pln.bixelWidth/2/0.1) = 1;

% 2D convolution of Fluence and Kernels in fourier domain
convMx1 = real(fftshift(ifft2(fft2(F).*fft2(kernel1Mx))));
convMx2 = real(fftshift(ifft2(fft2(F).*fft2(kernel2Mx))));
convMx3 = real(fftshift(ifft2(fft2(F).*fft2(kernel3Mx))));

% Creates an interpolant for kernes from vectors position X and Z
if exist('griddedInterpolant','class') % use griddedInterpoland class when available 
    Interp_kernel1 = griddedInterpolant(X',Z',convMx1','linear');
    Interp_kernel2 = griddedInterpolant(X',Z',convMx2','linear');
    Interp_kernel3 = griddedInterpolant(X',Z',convMx3','linear');
else
    Interp_kernel1 = @(x,y)interp2(X,Z,convMx1,x,y,'linear');
    Interp_kernel2 = @(x,y)interp2(X,Z,convMx2,x,y,'linear');
    Interp_kernel3 = @(x,y)interp2(X,Z,convMx3,x,y,'linear');
end

% generate meshgrid with CT position [mm]
[X_geo,Y_geo,Z_geo] = meshgrid(ct.resolution(1)*(1:size(ct.cube,2)),...
    ct.resolution(2)*(1:size(ct.cube,1)),ct.resolution(3)*(1:size(ct.cube,3)));

% take only voxels inside patient
X_geo = X_geo(V);
Y_geo = Y_geo(V);
Z_geo = Z_geo(V);

% define source position for beam eye view.
sourcePoint_bev = [0 -pln.SAD 0];

counter = 0;

fprintf('matRad: Photon dose calculation... ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams; % loop over all beams
    
    % gantry and couch roation matrices according to IEC 61217 standard
    % instead of moving the beam around the patient, we perform an inverse
    % rotation of the patient, i.e. we consider a beam's eye view
    % coordinate system
    
    % rotation around Z axis (gantry)
    rotMx_XY = [cosd(pln.gantryAngles(i)) -sind(pln.gantryAngles(i)) 0;
                sind(pln.gantryAngles(i))  cosd(pln.gantryAngles(i)) 0;
                                        0                          0 1];
    
    % rotation around Y axis (couch)
    rotMx_XZ = [ cosd(pln.couchAngles(i)) 0 sind(pln.couchAngles(i));
                                        0 1                         0;
                -sind(pln.couchAngles(i)) 0 cosd(pln.couchAngles(i))];
    
    % rotate target coordinates around Y axis and then around Z axis
    % i.e. 1st couch, 2nd gantry; matrix multiplication not cummutative
    rot_coords = coords_inside*rotMx_XZ*rotMx_XY;
    
    for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
        
        counter = counter + 1;

        % Display progress
        matRad_progress(counter,dij.totalNumOfBixels);
        waitbar(counter/dij.totalNumOfBixels);
        % remember beam and bixel number
        dij.beamNum(counter)  = i;
        dij.rayNum(counter)   = j;
        dij.bixelNum(counter) = j;
        
        % Ray tracing for beam i and bixel j
        [ix,radDepths,geoDists,latDistsX,latDistsZ] = matRad_calcRadGeoDists(ct.cube,V,...
            pln.isoCenter,rot_coords,ct.resolution,stf(i).sourcePoint,...
            stf(i).ray(j).targetPoint,sourcePoint_bev,...
            stf(i).ray(j).targetPoint_bev,X_geo,Y_geo,Z_geo,lateralCutoff,visBool);
        
        % calculate photon dose for beam i and bixel j
        bixelDose = matRad_calcPhotonDoseBixel(pln.SAD,m,betas, ...
                                               Interp_kernel1,...
                                               Interp_kernel2,...
                                               Interp_kernel3,...
                                               radDepths,...
                                               geoDists,...
                                               latDistsX,...
                                               latDistsZ);
       
        % Save dose for every bixel in cell array
        doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(V(ix),1,bixelDose,numel(ct.cube),1);
                
        % save computation time and memory by sequentially filling the 
        % sparse matrix dose.dij from the cell array
        if mod(counter,numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
            dij.physicalDose(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [doseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
        end
        
    end
end

close(figureWait);