function dij = matRad_calcPhotonDose(ct,stf,pln,cst,multScen)
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
%   multScen:   matRad multiple scnerio struct
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

% set up arrays for book keeping
dij.bixelNum = NaN*ones(dij.totalNumOfRays,1);
dij.rayNum   = NaN*ones(dij.totalNumOfRays,1);
dij.beamNum  = NaN*ones(dij.totalNumOfRays,1);

% Allocate space for dij.physicalDose sparse matrix
for CtScen = 1:multScen.numOfCtScen
    for ShiftScen = 1:multScen.numOfShiftScen
        for RangeShiftScen = 1:multScen.numOfRangeShiftScen  
            
            if multScen.ScenCombMask(CtScen,ShiftScen,RangeShiftScen)
                dij.physicalDose{CtScen,ShiftScen,RangeShiftScen} = spalloc(prod(ct.cubeDim),dij.totalNumOfBixels,1);
            end
            
        end
    end
end

% Allocate memory for dose_temp cell array
numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
doseTmpContainer = cell(numOfBixelsContainer,multScen.numOfCtScen,multScen.numOfShiftScen,multScen.numOfRangeShiftScen);

% take only voxels inside patient
V = [cst{:,4}];
V = unique(vertcat(V{:}));

% Convert CT subscripts to linear indices.
[yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(ct.cubeDim,V);

% set lateral cutoff value
lateralCutoff = 30; % [mm]

% toggle custom primary fluence on/off. if 0 we assume a homogeneous
% primary fluence, if 1 we use measured radially symmetric data
useCustomPrimFluenceBool = 0;

%% kernel convolution
% load polynomial fits for kernels ppKernel1, ppKernel2, ppKernel3
fileName = [pln.radiationMode '_' pln.machine];
try
   load(fileName);
catch
   error(['Could not find the following machine file: ' fileName ]); 
end

% Make a 2D grid extending +/-100mm with 0.1 mm resolution
convLimits = 100; % [mm]
convResolution = .5; % [mm]
[X,Z] = meshgrid(-convLimits:convResolution:convLimits);
                          
% Evaluate piecewise polynomial kernels
kernel1Mx = ppval(machine.data.ppKernel1,sqrt(X.^2+Z.^2));
kernel2Mx = ppval(machine.data.ppKernel2,sqrt(X.^2+Z.^2));
kernel3Mx = ppval(machine.data.ppKernel3,sqrt(X.^2+Z.^2));

% Create zero matrix for the Fluence
F = zeros(size(X));

% set bixel opening to one
F(abs(X)<=pln.bixelWidth/2 & abs(Z)<=pln.bixelWidth/2) = 1;

% gaussian convolution of field to model penumbra
sigmaGauss = 2.1/convResolution; % [mm] / see diploma thesis siggel 4.1.2
gaussFilter =  convResolution^2/(2*pi*sigmaGauss^2) * exp( -(X.^2+Z.^2)/(2*sigmaGauss^2) );
F = real(fftshift(ifft2(fft2( ifftshift(F) ).*fft2( ifftshift(gaussFilter) ))));

if ~useCustomPrimFluenceBool % pre-compute konvolution matrices for idealized homogeneous primary fluence
    
    % Display console message.
    fprintf('matRad: Uniform primary photon fluence -> pre-compute kernel convolution... \n');    

    % 2D convolution of Fluence and Kernels in fourier domain
    convMx1 = real(fftshift(ifft2(fft2( ifftshift(F) ).*fft2( ifftshift(kernel1Mx) ) )));
    convMx2 = real(fftshift(ifft2(fft2( ifftshift(F) ).*fft2( ifftshift(kernel2Mx) ) )));
    convMx3 = real(fftshift(ifft2(fft2( ifftshift(F) ).*fft2( ifftshift(kernel3Mx) ) )));

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

end

counter = 0;

fprintf('matRad: Photon dose calculation...\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams; % loop over all beams
    
    fprintf(['Beam ' num2str(i) ' of ' num2str(dij.numOfBeams) ': \n']);

    bixelsPerBeam = 0;

    % convert voxel indices to real coordinates using iso center of beam i
    for ShiftScen = 1:multScen.numOfShiftScen
        xCoordsV = xCoordsV_vox(:)*ct.resolution.x - (stf(i).isoCenter(1) + multScen.shifts(1,ShiftScen));
        yCoordsV = yCoordsV_vox(:)*ct.resolution.y - (stf(i).isoCenter(2) + multScen.shifts(2,ShiftScen));
        zCoordsV = zCoordsV_vox(:)*ct.resolution.z - (stf(i).isoCenter(3) + multScen.shifts(3,ShiftScen));
        coordsV{ShiftScen}  = [xCoordsV yCoordsV zCoordsV];
    end
    
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
    for ShiftScen = 1:multScen.numOfShiftScen
        rot_coordsV{ShiftScen} = coordsV{ShiftScen}*inv_rotMx_XZ_T*inv_rotMx_XY_T;
    
        rot_coordsV{ShiftScen}(:,1) = rot_coordsV{ShiftScen}(:,1)-stf(i).sourcePoint_bev(1);
        rot_coordsV{ShiftScen}(:,2) = rot_coordsV{ShiftScen}(:,2)-stf(i).sourcePoint_bev(2);
        rot_coordsV{ShiftScen}(:,3) = rot_coordsV{ShiftScen}(:,3)-stf(i).sourcePoint_bev(3);
    end
    
    % ray tracing
    fprintf(['matRad: calculate radiological depth cube...']);
    [radDepthCube,geoDistCube] = matRad_rayTracing(stf(i),ct,V,lateralCutoff,multScen);
    fprintf('done \n');
    
    % construct binary mask where ray tracing results are available
    for ShiftScen = 1:multScen.numOfShiftScen
        radDepthMask{ShiftScen} = ~isnan(radDepthCube{1,ShiftScen});
    end
    %radDepthIx = true(ct.cubeDim);                              % für ctScen überflüssig
    %for CtScen = 1:multScen.numOfCtScen                         % für ctScen überflüssig
    %    radDepthIx = radDepthIx .* isnan(radDepthCube{CtScen}); % für ctScen überflüssig
    %end                                                         % für ctScen überflüssig
    %radDepthIx = ~radDepthIx;                                   % für ctScen überflüssig

    for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
        
        counter = counter + 1;
        bixelsPerBeam = bixelsPerBeam + 1;
        
        if useCustomPrimFluenceBool % use custom primary fluence if specifried
            
            r     = sqrt( (X-stf(i).ray(j).rayPos(1)).^2 + (Z-stf(i).ray(j).rayPos(3)).^2 );
            Psi   = interp1(primaryFluence(:,1),primaryFluence(:,2),r);
            FxPsi = F .* Psi;
        
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

        end

        % Display progress
        matRad_progress(bixelsPerBeam,stf(i).totalNumOfBixels);
        % update waitbar only 100 times
        if mod(counter,round(dij.totalNumOfBixels/100)) == 0 && ishandle(figureWait)
            waitbar(counter/dij.totalNumOfBixels);
        end
        
        % remember beam and bixel number
        dij.beamNum(counter)  = i;
        dij.rayNum(counter)   = j;
        dij.bixelNum(counter) = j;
        
        
        % Ray tracing for beam i and bixel j
        for ShiftScen = 1:multScen.numOfShiftScen
            [ix{ShiftScen},~,latDistsX{ShiftScen},latDistsZ{ShiftScen}] = matRad_calcGeoDists(rot_coordsV{ShiftScen}, ...
                                                                          stf(i).sourcePoint_bev, ...
                                                                          stf(i).ray(j).targetPoint_bev, ...
                                                                          geoDistCube{ShiftScen}(V), ...
                                                                          machine.meta.SAD, ...
                                                                          radDepthMask{ShiftScen}(V), ...
                                                                          lateralCutoff);
        end
                                                       
        for CtScen = 1:multScen.numOfCtScen
            for ShiftScen = 1:multScen.numOfShiftScen
                for RangeShiftScen = 1:multScen.numOfRangeShiftScen  
            
                    if multScen.ScenCombMask(CtScen,ShiftScen,RangeShiftScen)
                        
                        % manipulate radDepthCube for range scenarios
                        radDepths = radDepthCube{CtScen,ShiftScen}(V(ix{ShiftScen}));                                         

                        if multScen.relRangeShifts(RangeShiftScen) ~= 0 || multScen.absRangeShifts(RangeShiftScen) ~= 0
                            radDepths = radDepths +...                                                                                % original cube
                                        radDepthCube{CtScen,ShiftScen}(V(ix{ShiftScen}))*multScen.relRangeShifts(RangeShiftScen) +... % rel range shift
                                        multScen.absRangeShifts(RangeShiftScen);                                                      % absolute range shift
                            radDepths(radDepths < 0) = 0;  
                        end                    
                                              
                        % calculate photon dose for beam i and bixel j
                        bixelDose = matRad_calcPhotonDoseBixel(machine.meta.SAD,machine.data.m,...
                                                               machine.data.betas, ...
                                                               Interp_kernel1,...
                                                               Interp_kernel2,...
                                                               Interp_kernel3,...
                                                               manipulatedRadDepthCube,...
                                                               geoDistCube{ShiftScen}(V(ix{ShiftScen})),...
                                                               latDistsX{ShiftScen},...
                                                               latDistsZ{ShiftScen});

                        % Save dose for every bixel in cell array
                        doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,CtScen,ShiftScen,RangeShiftScen} = sparse(V(ix{ShiftScen}),1,bixelDose,prod(ct.cubeDim),1);
                    end
                    
                end
            end
        end
                
        % save computation time and memory by sequentially filling the 
        % sparse matrix dose.dij from the cell array
        if mod(counter,numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
            for CtScen = 1:multScen.numOfCtScen
                for ShiftScen = 1:multScen.numOfShiftScen
                    for RangeShiftScen = 1:multScen.numOfRangeShiftScen
                        
                        if multScen.ScenCombMask(CtScen,ShiftScen,RangeShiftScen)
                            dij.physicalDose{CtScen,ShiftScen,RangeShiftScen}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [doseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,CtScen,ShiftScen,RangeShiftScen}];
                        end
                        
                    end
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


