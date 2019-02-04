function dij = matRad_calcParticleDoseMC(ct,stf,pln,cst,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad MCsqaure monte carlo photon dose calculation wrapper
%
% call
%   dij = matRad_calcParticleDoseMc(ct,stf,pln,cst,calcDoseDirect)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   calcDoseDirect:             binary switch to enable forward dose
%                               calcualtion
% output
%   dij:                        matRad dij struct
%
% References
%
%   https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.4943377
%   http://www.openmcsquare.org/
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if valid machine
if ~strcmp(pln.radiationMode,'protons') || ~strcmp(pln.machine,'generic_MCsquare')
    error('wrong radiation modality and/or machine.');    
end

%% check if binaries are available
if ispc
    if exist('MCSquare_windows.exe') ~= 2
        error('Could not find MCsquare binary.\n');
    end
elseif isunix
    error('xxxxxxxxxxxxxx\n');
elseif ismac
    error('MCsquare binaries not available for mac OS.\n');
end

if nargin < 5
    calcDoseDirect = false;
end

% set and change to MCsquare binary folder
currFolder = pwd;
fullfilename = mfilename('fullpath');
MCsquareFolder = [fullfilename(1:find(fullfilename==filesep,1,'last')) 'submodules' filesep 'MCsquare'];

% cd to MCsquare folder (necessary for binary)
cd(MCsquareFolder);

% to guarantee downwards compatibility with data that does not have
% ct.x/y/z
if ~any(isfield(ct,{'x','y','z'}))
    ct.x = ct.resolution.x*[1:ct.cubeDim(1)];
    ct.y = ct.resolution.y*[1:ct.cubeDim(2)];
    ct.z = ct.resolution.z*[1:ct.cubeDim(3)];
end

% set grids
if ~isfield(pln,'propDoseCalc') || ...
   ~isfield(pln.propDoseCalc,'doseGrid') || ...
   ~isfield(pln.propDoseCalc.doseGrid,'resolution')
    % default values
    dij.doseGrid.resolution.x = 2.5; % [mm]
    dij.doseGrid.resolution.y = 2.5; % [mm]
    dij.doseGrid.resolution.z = 2.5;   % [mm]
else
    % take values from pln strcut
    dij.doseGrid.resolution.x = pln.propDoseCalc.doseGrid.resolution.x;
    dij.doseGrid.resolution.y = pln.propDoseCalc.doseGrid.resolution.y;
    dij.doseGrid.resolution.z = pln.propDoseCalc.doseGrid.resolution.z;
end

% check if using isotropic dose grid resolution in x and y direction
if dij.doseGrid.resolution.x ~= dij.doseGrid.resolution.y
    error('Anisotropic resolution in x and y direction');
end

dij.doseGrid.x = ct.x(1):dij.doseGrid.resolution.x:ct.x(end);
dij.doseGrid.y = ct.y(1):dij.doseGrid.resolution.y:ct.y(end);
dij.doseGrid.z = ct.z(1):dij.doseGrid.resolution.z:ct.z(end);

dij.doseGrid.dimensions  = [numel(dij.doseGrid.y) numel(dij.doseGrid.x) numel(dij.doseGrid.z)];
dij.doseGrid.numOfVoxels = prod(dij.doseGrid.dimensions);

dij.ctGrid.resolution.x = ct.resolution.x;
dij.ctGrid.resolution.y = ct.resolution.y;
dij.ctGrid.resolution.z = ct.resolution.z;

dij.ctGrid.x = ct.x;
dij.ctGrid.y = ct.y;
dij.ctGrid.z = ct.z;

dij.ctGrid.dimensions  = [numel(dij.ctGrid.y) numel(dij.ctGrid.x) numel(dij.ctGrid.z)];
dij.ctGrid.numOfVoxels = prod(dij.ctGrid.dimensions);

% meta information for dij
dij.numOfBeams         = pln.propStf.numOfBeams;
dij.numOfScenarios     = 1;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
if calcDoseDirect
    dij.totalNumOfBixels = 1;
else
    dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
end
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);

% book keeping
dij.bixelNum          = NaN*ones(dij.totalNumOfBixels,1);
dij.rayNum            = NaN*ones(dij.totalNumOfBixels,1);
dij.beamNum           = NaN*ones(dij.totalNumOfBixels,1);
dij.MCsquareCalcOrder = NaN*ones(dij.totalNumOfBixels,1);

% take only voxels inside patient
VctGrid = [cst{:,4}];
VctGrid = unique(vertcat(VctGrid{:}));

% receive linear indices and grid locations from the dose grid
tmpCube    = zeros(ct.cubeDim);
tmpCube(VctGrid) = 1;
% interpolate cube
VdoseGrid = find(interp3(dij.ctGrid.x,  dij.ctGrid.y,   dij.ctGrid.z,tmpCube, ...
                         dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest'));

% ignore densities outside of contours
eraseCtDensMask = ones(dij.ctGrid.numOfVoxels,1);
eraseCtDensMask(VctGrid) = 0;
for i = 1:ct.numOfCtScen
    ct.cubeHU{i}(eraseCtDensMask == 1) = -1024;
end

% downsample ct
for s = 1:dij.numOfScenarios
    HUcube{s} =  interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cubeHU{s}, ...
                         dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
end

% what are you doing here, Lucas?
nbThreads = 4;

% set number of particles simulated per pencil beam
if calcDoseDirect
    nCasePerBixel = 1000000;
else
    nCasePerBixel = 100000;
end

% set relative dose cutoff for storage in dose influence matrix
relDoseCutoff = 10^(-4);
% set absolute calibration factor
% convert from eV/g/primary to Gy 1e6 primaries
absCalibrationFactorMC2 = 1.602176e-19 * 1.0e+9;

% Allocate space for dij.physicalDose sparse matrix
for i = 1:dij.numOfScenarios
    dij.physicalDose{i} = spalloc(dij.doseGrid.numOfVoxels,dij.totalNumOfBixels,1);
end

if ~strcmp(pln.radiationMode,'protons')
    errordlg('MCsquare is only supported for protons');
end

if isequal(pln.propOpt.bioOptimization,'const_RBExD')
            dij.RBE = 1.1;
            fprintf(['matRad: Using a constant RBE of 1.1 \n']);
end

% MCsquare settings
MCsquareConfigFile = 'MCsquareConfig.txt';

MCsquareConfig = MatRad_MCsquareConfig;

MCsquareConfig.BDL_Plan_File = 'currBixels.txt';
MCsquareConfig.CT_File       = 'MC2patientCT.mhd';
MCsquareConfig.Num_Threads   = nbThreads;
MCsquareConfig.RNG_Seed      = 1234;
MCsquareConfig.Num_Primaries = nCasePerBixel;

% turn simulation of individual beamlets
MCsquareConfig.Beamlet_Mode = ~calcDoseDirect;
% turn of writing of full dose cube
MCsquareConfig.Dose_MHD_Output = calcDoseDirect;
% turn on sparse output
MCsquareConfig.Dose_Sparse_Output = ~calcDoseDirect;
% set threshold of sparse matrix generation
MCsquareConfig.Dose_Sparse_Threshold = relDoseCutoff;

% write patient data
MCsquareBinCubeResolution = [dij.doseGrid.resolution.x ...
                             dij.doseGrid.resolution.y ...
                             dij.doseGrid.resolution.z];   
matRad_writeMhd(HUcube{1},MCsquareBinCubeResolution,MCsquareConfig.CT_File);

% prepare steering for MCsquare and sort stf by energy
isoCenterOffset = [dij.doseGrid.resolution.x/2 dij.doseGrid.resolution.y/2 dij.doseGrid.resolution.z/2] ...
                - [dij.ctGrid.resolution.x   dij.ctGrid.resolution.y   dij.ctGrid.resolution.z];

counter = 0;             
for i = 1:length(stf)
    stfMCsquare(i).gantryAngle = mod(180-stf(i).gantryAngle,360);
    stfMCsquare(i).couchAngle  = stf(i).couchAngle;
    stfMCsquare(i).isoCenter   = stf(i).isoCenter + isoCenterOffset;
    stfMCsquare(i).energies    = unique([stf(i).ray.energy]);
    
    % allocate empty target point container
    for j = 1:numel(stfMCsquare(i).energies)
        stfMCsquare(i).energyLayer(j).targetPoints   = [];
        stfMCsquare(i).energyLayer(j).numOfPrimaries = [];
        stfMCsquare(i).energyLayer(j).rayNum         = [];
        stfMCsquare(i).energyLayer(j).bixelNum       = [];
    end
    
    for j = 1:stf(i).numOfRays
        for k = 1:stf(i).numOfBixelsPerRay(j)
            counter = counter + 1;
            dij.beamNum(counter)  = i;
            dij.rayNum(counter)   = j;
            dij.bixelNum(counter) = k;
        end
        
        for k = 1:numel(stfMCsquare(i).energies)
            if any(stf(i).ray(j).energy == stfMCsquare(i).energies(k))
                stfMCsquare(i).energyLayer(k).rayNum   = [stfMCsquare(i).energyLayer(k).rayNum j];
                stfMCsquare(i).energyLayer(k).bixelNum = [stfMCsquare(i).energyLayer(k).bixelNum ...
                    find(stf(i).ray(j).energy == stfMCsquare(i).energies(k))];
                stfMCsquare(i).energyLayer(k).targetPoints = [stfMCsquare(i).energyLayer(k).targetPoints; ...
                                        -stf(i).ray(j).rayPos_bev(1) stf(i).ray(j).rayPos_bev(3)];
                if calcDoseDirect
                    stfMCsquare(i).energyLayer(k).numOfPrimaries = [stfMCsquare(i).energyLayer(k).numOfPrimaries ...
                                         round(stf(i).ray(j).weight(stf(i).ray(j).energy == stfMCsquare(i).energies(k))*MCsquareConfig.Num_Primaries)];
                else
                    stfMCsquare(i).energyLayer(k).numOfPrimaries = [stfMCsquare(i).energyLayer(k).numOfPrimaries ...
                        MCsquareConfig.Num_Primaries];
                end
            end
        end
               
    end
    
end

% remember order
counterMCsquare = 0;
MCsquareOrder = NaN * ones(dij.totalNumOfBixels,1);
for i = 1:length(stf)
    for j = 1:numel(stfMCsquare(i).energies)
        for k = 1:numel(stfMCsquare(i).energyLayer(j).numOfPrimaries)
            counterMCsquare = counterMCsquare + 1;
            ix = find(i                                         == dij.beamNum & ...
                      stfMCsquare(i).energyLayer(j).rayNum(k)   == dij.rayNum & ...
                      stfMCsquare(i).energyLayer(j).bixelNum(k) == dij.bixelNum);
        
            MCsquareOrder(ix) = counterMCsquare;
        end
    end
end

if any(isnan(MCsquareOrder))
    error('Wrong order')
end

%% MC computation and dij filling
matRad_writeMCsquareinputAllFiles(MCsquareConfigFile,MCsquareConfig,stfMCsquare);

% run MCsquare
[status,cmdout] = system(['MCSquare_windows.exe ' MCsquareConfigFile],'-echo');

mask = false(dij.doseGrid.numOfVoxels,1);
mask(VdoseGrid) = true;

% read sparse matrix
if ~calcDoseDirect
    dij.physicalDose{1} = absCalibrationFactorMC2 * matRad_sparseBeamletsReaderMCsquare ( ...
                    [MCsquareConfig.Output_Directory filesep 'Sparse_Dose.bin'], ...
                    dij.doseGrid.dimensions, ...
                    dij.totalNumOfBixels, ...
                    mask);
else
    cube = matRad_readMhd(MCsquareConfig.Output_Directory,'Dose.mhd');
    dij.physicalDose{1} = sparse(VdoseGrid,ones(numel(VdoseGrid),1), ...
                                 absCalibrationFactorMC2 * cube(VdoseGrid), ...
                                 dij.doseGrid.numOfVoxels,1);
end

% reorder influence matrix to comply with matRad default ordering
if MCsquareConfig.Beamlet_Mode
    dij.physicalDose{1} = dij.physicalDose{1}(:,MCsquareOrder);            
end        

fprintf('matRad: done!\n');

try
    % wait 0.1s for closing all waitbars
    allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar');
    delete(allWaitBarFigures);
    pause(0.1);
catch
end

%% clear all data
delete([MCsquareConfig.CT_File(1:end-4) '.*']);
delete('currBixels.txt');
delete('MCsquareConfig.txt');
eval(['rmdir ' MCsquareConfig.Output_Directory ' s']);

% cd back
cd(currFolder);

end