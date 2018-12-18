function dij = matRad_calcPhotonDose(ct,stf,pln,cst,param)
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
%   param:          (optional) structure defining additional parameter
%                   param.calcDoseDirect boolian switch to bypass dose influence matrix
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
% Copyright 2017 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if exist('param','var')
   if ~isfield(param,'logLevel')
      param.logLevel = 1;
   end
   % default: dose influence matrix computation
   if ~isfield(param,'calcDoseDirect')
      param.calcDoseDirect = false;
   end
   
else
   param.calcDoseDirect = false;
   param.subIx          = [];
   param.logLevel       = 1;
end

% set consistent random seed (enables reproducibility)
[env, ~] = matRad_getEnvironment();

switch env
   case 'MATLAB'
      rng(0);
   case 'OCTAVE'
      rand('seed',0)
end

if param.logLevel == 1
   % initialize waitbar
   figureWait = waitbar(0,'calculate dose influence matrix for photons...');
   % show busy state
   set(figureWait,'pointer','watch');
end

% to guarantee downwards compatibility with data that does not have
% ct.x/y/z
if ~any(isfield(ct,{'x','y','z'}))
   ct.x = ct.resolution.x*[0:ct.cubeDim(1)-1]-ct.resolution.x/2;
   ct.y = ct.resolution.y*[0:ct.cubeDim(2)-1]-ct.resolution.y/2;
   ct.z = ct.resolution.z*[0:ct.cubeDim(3)-1]-ct.resolution.z/2;
end

% set grids
if ~isfield(pln,'propDoseCalc') || ...
      ~isfield(pln.propDoseCalc,'doseGrid') || ...
      ~isfield(pln.propDoseCalc.doseGrid,'resolution')
   % default values
   dij.doseGrid.resolution.x = 2.5; % [mm]
   dij.doseGrid.resolution.y = 2.5; % [mm]
   dij.doseGrid.resolution.z = 3;   % [mm]
else
   % take values from pln strcut
   dij.doseGrid.resolution.x = pln.propDoseCalc.doseGrid.resolution.x;
   dij.doseGrid.resolution.y = pln.propDoseCalc.doseGrid.resolution.y;
   dij.doseGrid.resolution.z = pln.propDoseCalc.doseGrid.resolution.z;
end

dij.doseGrid.x = ct.x(1):dij.doseGrid.resolution.x:ct.x(end);
dij.doseGrid.y = ct.y(1):dij.doseGrid.resolution.y:ct.y(end);
dij.doseGrid.z = ct.z(1):dij.doseGrid.resolution.z:ct.z(end);

dij.doseGrid.dimensions  = [numel(dij.doseGrid.x) numel(dij.doseGrid.y) numel(dij.doseGrid.z)];
dij.doseGrid.numOfVoxels = prod(dij.doseGrid.dimensions);

dij.ctGrid.resolution.x = ct.resolution.x;
dij.ctGrid.resolution.y = ct.resolution.y;
dij.ctGrid.resolution.z = ct.resolution.z;

dij.ctGrid.x = ct.x;
dij.ctGrid.y = ct.y;
dij.ctGrid.z = ct.z;

dij.ctGrid.dimensions  = [numel(dij.ctGrid.x) numel(dij.ctGrid.y) numel(dij.ctGrid.z)];
dij.ctGrid.numOfVoxels = prod(dij.ctGrid.dimensions);

% calculate rED or rSP from HU
ct = matRad_calcWaterEqD(ct, pln, param);

% meta information for dij
dij.numOfBeams         = pln.propStf.numOfBeams;
dij.numOfScenarios     = pln.multScen.numOfCtScen;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);

% check if full dose influence data is required
if param.calcDoseDirect
   numOfColumnsDij      = length(stf);
   numOfBixelsContainer = 1;
else
   numOfColumnsDij      = dij.totalNumOfBixels;
   numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
end

% set up arrays for book keeping
dij.bixelNum = NaN*ones(numOfColumnsDij,1);
dij.rayNum   = NaN*ones(numOfColumnsDij,1);
dij.beamNum  = NaN*ones(numOfColumnsDij,1);

% Allocate space for dij.physicalDose sparse matrix
for CtScen = 1:pln.multScen.numOfCtScen
   for ShiftScen = 1:pln.multScen.totNumShiftScen
      for RangeShiftScen = 1:pln.multScen.totNumRangeScen
         
         if pln.multScen.scenMask(CtScen,ShiftScen,RangeShiftScen)
            dij.physicalDose{CtScen,ShiftScen,RangeShiftScen} = spalloc(dij.doseGrid.numOfVoxels,numOfColumnsDij,1);
         end
         
      end
   end
end

doseTmpContainer = cell(numOfBixelsContainer,pln.multScen.numOfCtScen,pln.multScen.totNumShiftScen,pln.multScen.totNumRangeScen);

% Only take voxels inside patient at full ct resolution
if isfield(param, 'subIx') && ~isempty(param.subIx) && param.calcDoseDirect
   VctGrid = param.subIx;
else
   VctGrid = [cst{:,4}];
   VctGrid = unique(vertcat(VctGrid{:}));
end

% ignore densities outside of contours
eraseCtDensMask = ones(prod(ct.cubeDim),1);
eraseCtDensMask(VctGrid) = 0;
for i = 1:ct.numOfCtScen
   ct.cube{i}(eraseCtDensMask == 1) = 0;
end

% ser overlap prioriites
cst = matRad_setOverlapPriorities(cst);

% resizing cst to dose cube resolution
cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,...
   dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);

% Convert CT subscripts to linear indices.
[yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(ct.cubeDim,VctGrid);

% receive linear indices and grid locations from the dose grid
tmpCube    = zeros(ct.cubeDim);
tmpCube(VctGrid) = 1;
% interpolate cube
VdoseGrid = find(interp3(dij.ctGrid.y,  dij.ctGrid.x,   dij.ctGrid.z,tmpCube, ...
   dij.doseGrid.y,dij.doseGrid.x',dij.doseGrid.z)>0.5);

% Convert CT subscripts to coarse linear indices.
[yCoordsV_voxDoseGrid, xCoordsV_voxDoseGrid, zCoordsV_voxDoseGrid] = ind2sub(dij.doseGrid.dimensions,VdoseGrid);

% set lateral cutoff value
lateralCutoff = 50; % [mm]

% toggle custom primary fluence on/off. if 0 we assume a homogeneous
% primary fluence, if 1 we use measured radially symmetric data
useCustomPrimFluenceBool = 0;

% 0 if field calc is bixel based, 1 if dose calc is field based
isFieldBasedDoseCalc = strcmp(num2str(pln.propStf.bixelWidth),'field');

%% kernel convolution
% prepare data for convolution to reduce calculation time
fileName = [pln.radiationMode '_' pln.machine];
try
   load([fileparts(mfilename('fullpath')) filesep fileName]);
catch
   matRad_dispToConsole(['Could not find the following machine file: ' fileName ],param,'error');
end

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

% book keeping - this is necessary since pln is not used in optimization or
% matRad_calcCubes
if strcmp(pln.bioParam.model,'constRBE')
   dij.RBE = pln.bioParam.RBE;
end

for ShiftScen = 1:pln.multScen.totNumShiftScen
   
   % manipulate isocenter
   for k = 1:numel(stf)
      stf(k).isoCenter = stf(k).isoCenter + pln.multScen.isoShift(ShiftScen,:);
   end
   
   counter = 0;
   % compute SSDs
   stf = matRad_computeSSD(stf,ct);
   
   matRad_dispToConsole(['shift scenario ' num2str(ShiftScen) ' of ' num2str(pln.multScen.totNumShiftScen) ': \n'],param,'info');
   matRad_dispToConsole('matRad: photon dose calculation...\n',param,'info');
   
   for i = 1:numel(stf) % loop over all beams
      
      matRad_dispToConsole(['Beam ' num2str(i) ' of ' num2str(dij.numOfBeams) ': \n'],param,'info');
      
      % remember beam and bixel number
      if param.calcDoseDirect
         dij.beamNum(i)    = i;
         dij.rayNum(i)     = i;
         dij.bixelNum(i)   = i;
      end
      
      bixelsPerBeam = 0;
      
      % convert voxel indices to real coordinates using iso center of beam i
      xCoordsV       = xCoordsV_vox(:)*ct.resolution.x-stf(i).isoCenter(1);
      yCoordsV       = yCoordsV_vox(:)*ct.resolution.y-stf(i).isoCenter(2);
      zCoordsV       = zCoordsV_vox(:)*ct.resolution.z-stf(i).isoCenter(3);
      coordsV        = [xCoordsV yCoordsV zCoordsV];
      
      xCoordsVdoseGrid = xCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.x-stf(i).isoCenter(1);
      yCoordsVdoseGrid = yCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.y-stf(i).isoCenter(2);
      zCoordsVdoseGrid = zCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.z-stf(i).isoCenter(3);
      coordsVdoseGrid  = [xCoordsVdoseGrid yCoordsVdoseGrid zCoordsVdoseGrid];
      
      % Get Rotation Matrix
      % Do not transpose matrix since we usage of row vectors &
      % transformation of the coordinate system need double transpose
      
      rotMat_system_T = matRad_getRotationMatrix(pln.propStf.gantryAngles(i),pln.propStf.couchAngles(i));
      
      
      % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
      rot_coordsV         = coordsV*rotMat_system_T;
      rot_coordsVdoseGrid = coordsVdoseGrid*rotMat_system_T;
      
      rot_coordsV(:,1) = rot_coordsV(:,1)-stf(i).sourcePoint_bev(1);
      rot_coordsV(:,2) = rot_coordsV(:,2)-stf(i).sourcePoint_bev(2);
      rot_coordsV(:,3) = rot_coordsV(:,3)-stf(i).sourcePoint_bev(3);
      
      rot_coordsVdoseGrid(:,1) = rot_coordsVdoseGrid(:,1)-stf(i).sourcePoint_bev(1);
      rot_coordsVdoseGrid(:,2) = rot_coordsVdoseGrid(:,2)-stf(i).sourcePoint_bev(2);
      rot_coordsVdoseGrid(:,3) = rot_coordsVdoseGrid(:,3)-stf(i).sourcePoint_bev(3);
      
      % calculate geometric distances
      geoDistVdoseGrid{1}= sqrt(sum(rot_coordsVdoseGrid.^2,2));
      
      % ray tracing
      matRad_dispToConsole('matRad: calculate radiological depth cube...',param,'info');
      radDepthVctGrid = matRad_rayTracing(stf(i),ct,VctGrid,rot_coordsV,effectiveLateralCutoff);
      matRad_dispToConsole('done \n',param,'info');
      
      % interpolate radiological depth cube to dose grid resolution
      radDepthVdoseGrid = matRad_interpRadDepth...
         (ct,ct.numOfCtScen,VctGrid,VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,radDepthVctGrid);
      
      % limit rotated coordinates to positions where ray tracing is availabe
      rot_coordsVdoseGrid = rot_coordsVdoseGrid(~isnan(radDepthVdoseGrid{1}),:);
      
      % get index of central ray or closest to the central ray
      [~,center] = min(sum(reshape([stf(i).ray.rayPos_bev],3,[]).^2));
      
      % get correct kernel for given SSD at central ray (nearest neighbor approximation)
      [~,currSSDIx] = min(abs([machine.data.kernel.SSD]-stf(i).ray(center).SSD));
      
      matRad_dispToConsole(['                   SSD = ' num2str(machine.data.kernel(currSSDIx).SSD) 'mm                 \n'],param,'info');
      
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
         
         counter       = counter + 1;
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
         
         if param.logLevel <= 2
            % Display progress and update text only 200 times
            if mod(bixelsPerBeam,max(1,round(stf(i).totalNumOfBixels/200))) == 0
               matRad_progress(bixelsPerBeam/max(1,round(stf(i).totalNumOfBixels/200)),...
                  floor(stf(i).totalNumOfBixels/max(1,round(stf(i).totalNumOfBixels/200))));
            end
            if param.logLevel == 1
               % update waitbar only 100 times
               if mod(counter,round(dij.totalNumOfBixels/100)) == 0 && ishandle(figureWait)
                  waitbar(counter/dij.totalNumOfBixels);
               end
            end
         end
         
         % remember beam and bixel number
         if ~param.calcDoseDirect
            dij.beamNum(counter)  = i;
            dij.rayNum(counter)   = j;
            dij.bixelNum(counter) = j;
         end
         
         for CtScen = 1:pln.multScen.numOfCtScen
            for RangeShiftScen = 1:pln.multScen.totNumRangeScen
               
               if pln.multScen.scenMask(CtScen,ShiftScen,RangeShiftScen)
                  
                  % Ray tracing for beam i and bixel j
                  [ix,rad_distancesSq,isoLatDistsX,isoLatDistsZ] = matRad_calcGeoDists(rot_coordsVdoseGrid, ...
                     stf(i).sourcePoint_bev, ...
                     stf(i).ray(j).targetPoint_bev, ...
                     machine.meta.SAD, ...
                     find(~isnan(radDepthVdoseGrid{CtScen})), ...
                     effectiveLateralCutoff);
                  
                  % empty bixels may happen during recalculation of error
                  % scenarios -> skip to next bixel
                  if isempty(ix)
                     continue;
                  end
                  
                  % manipulate radDepthCube for range scenarios
                  manipulatedRadDepthCube = radDepthVdoseGrid{CtScen}(ix);
                  
                  if pln.multScen.relRangeShift(RangeShiftScen) ~= 0 || pln.multScen.absRangeShift(RangeShiftScen) ~= 0
                     manipulatedRadDepthCube = manipulatedRadDepthCube +...                                                                                % original cube
                        radDepthVdoseGrid{CtScen}(ix)*pln.multScen.relRangeShift(RangeShiftScen) +... % rel range shift
                        pln.multScen.absRangeShift(RangeShiftScen);                                                      % absolute range shift
                     manipulatedRadDepthCube(manipulatedRadDepthCube < 0) = 0;
                  end
                  
                  % calculate photon dose for beam i and bixel j
                  bixelDose = matRad_calcPhotonDoseBixel(machine.meta.SAD,machine.data.m,...
                     machine.data.betas, ...
                     Interp_kernel1,...
                     Interp_kernel2,...
                     Interp_kernel3,...
                     manipulatedRadDepthCube,...
                     geoDistVdoseGrid{1}(ix),...
                     isoLatDistsX,...
                     isoLatDistsZ);
                  
                  % sample dose only for bixel based dose calculation
                  if ~isFieldBasedDoseCalc
                     r0   = 25;   % [mm] sample beyond the inner core
                     Type = 'radius';
                     [ix,bixelDose] = matRad_DijSampling(ix,bixelDose,manipulatedRadDepthCube,rad_distancesSq,Type,r0);
                  end
                  % Save dose for every bixel in cell array
                  doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,CtScen,ShiftScen,RangeShiftScen} = sparse(VdoseGrid(ix),1,bixelDose,dij.doseGrid.numOfVoxels,1);
               end
            end
         end
         
         
         % save computation time and memory by sequentially filling the
         % sparse matrix dose.dij from the cell array
         if mod(counter,numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
            for CtScen = 1:pln.multScen.numOfCtScen
               for RangeShiftScen = 1:pln.multScen.totNumRangeScen
                  
                  if pln.multScen.scenMask(CtScen,ShiftScen,RangeShiftScen)
                     if param.calcDoseDirect
                        if isfield(stf(1).ray(1),'weight')
                           % score physical dose
                           dij.physicalDose{CtScen,ShiftScen,RangeShiftScen}(:,i) = dij.physicalDose{CtScen,ShiftScen,RangeShiftScen}(:,i) + stf(i).ray(j).weight * doseTmpContainer{1,CtScen,ShiftScen,RangeShiftScen};
                        else
                           matRad_dispToConsole(['No weight available for beam ' num2str(i) ', ray ' num2str(j)],param,'error');
                        end
                     else
                        % fill entire dose influence matrix
                        dij.physicalDose{CtScen,ShiftScen,RangeShiftScen}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [doseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,CtScen,ShiftScen,RangeShiftScen}];
                     end
                  end
                  
               end
            end
         end
         
      end
   end
   
   % manipulate isocenter
   for k = 1:length(stf)
      stf(k).isoCenter = stf(k).isoCenter - pln.multScen.isoShift(ShiftScen,:);
   end
   
end

% remove dose influence for voxels outside of segmentations for every ct
% scenario
for i = 1:pln.multScen.numOfCtScen
   
   % generate index set to erase
   tmpIx = [];
   for j = 1:size(cst,1)
      tmpIx = unique([tmpIx; cst{j,4}{i}]);
   end
   ix = setdiff(1:dij.doseGrid.numOfVoxels,tmpIx);
   
   for j = 1:pln.multScen.totNumShiftScen
      for k = 1:pln.multScen.totNumRangeScen
         if pln.multScen.scenMask(i,j,k)
            
            dij.physicalDose{i,j,k}(ix,:)      = 0;
            
            if isfield(dij,'mLETDose')
               dij.mLETDose{i,j,k}(ix,:)      = 0;
            end
            
            if pln.bioParam.bioOpt
               dij.mAlphaDose{i,j,k}(ix,:)    = 0;
               dij.mSqrtBetaDose{i,j,k}(ix,:) = 0;
            end
            
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

