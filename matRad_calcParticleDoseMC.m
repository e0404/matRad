function dij = matRad_calcParticleDoseMC(ct,stf,pln,cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad MCsqaure monte carlo photon dose calculation wrapper
%
% call
%   dij = matRad_calcParticleDoseMc(ct,stf,pln,cst,visBool)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   visBool:                    binary switch to enable visualization
% output
%   dij:                        matRad dij struct
%
% References
%
%   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX there needs to go a lot of stuff!
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
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
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);

% set up arrays for book keeping
dij.bixelNum = NaN*ones(dij.totalNumOfBixels,1);
dij.rayNum   = NaN*ones(dij.totalNumOfBixels,1);
dij.beamNum  = NaN*ones(dij.totalNumOfBixels,1);

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
nCasePerBixel = 10000;

% set relative dose cutoff for storage in dose influence matrix
relDoseCutoff = 10^(-3);
% set absolute calibration factor
% convert from sum of dose in Gy for all histories to Gy per 1e6 primaries
absCalibrationFactorMC2 = 1.602176e-19 * nCasePerBixel * 1.0e+6;

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
MCsquareConfig.Beamlet_Mode = true;
% turn of writing of full dose cube
MCsquareConfig.Dose_MHD_Output = false;
% turn on sparse output
MCsquareConfig.Dose_Sparse_Output = true;
% set threshold of sparse matrix generation
MCsquareConfig.Dose_Sparse_Threshold = relDoseCutoff;

% write patient data
MCsquareBinCubeResolution = [dij.doseGrid.resolution.x ...
                             dij.doseGrid.resolution.y ...
                             dij.doseGrid.resolution.z];   
matRad_writeMhd(HUcube{1},MCsquareBinCubeResolution,MCsquareConfig.CT_File);

% prepare steering for MCsquare and sort stf by energy
isoCenterOffset = [dij.doseGrid.resolution.x*1.5 dij.doseGrid.resolution.y/2 dij.doseGrid.resolution.z/2] ...
                - [dij.ctGrid.resolution.x   dij.ctGrid.resolution.y   dij.ctGrid.resolution.z];
            
for i = 1:length(stf)
    stfMCsquare(i).gantryAngle = mod(180-stf(i).gantryAngle,360);
    stfMCsquare(i).couchAngle  = stf(i).couchAngle;
    stfMCsquare(i).isoCenter   = stf(i).isoCenter + isoCenterOffset;
    stfMCsquare(i).energies    = unique([stf(i).ray.energy]);
    
    % allocate empty target point container
    for j = 1:numel(stfMCsquare(i).energies)
        stfMCsquare(i).energyLayer(j).targetPoints = [];
    end
    
    for j = 1:stf(i).numOfRays
        for k = 1:numel(stfMCsquare(i).energies)
            if any(stf(i).ray(j).energy == stfMCsquare(i).energies(k))
                stfMCsquare(i).energyLayer(k).targetPoints = [stfMCsquare(i).energyLayer(k).targetPoints; ...
                                        -stf(i).ray(j).rayPos_bev(1) stf(i).ray(j).rayPos_bev(3)];
            end
        end
               
    end
    
end

%% MC computation and dij filling
matRad_writeMCsquareinputAllFiles(MCsquareConfigFile,MCsquareConfig,stfMCsquare);



% run MCsquare
[status,cmdout] = system(['MCSquare_windows.exe ' MCsquareConfigFile],'-echo');


binSparseFileHeader = Sparse_read_header([MCsquareConfig.Output_Directory filesep ...
                                            'Sparse_Dose.txt']);
                    
mask = false(prod(binSparseFileHeader.ImageSize),1);
mask(VdoseGrid) = true;

dij.physicalDose{1} = absCalibrationFactorMC2 * matRad_sparseBeamletsReaderMSsquare ( ...
                [MCsquareConfig.Output_Directory filesep 'Sparse_Dose.bin'], ...
                binSparseFileHeader.ImageSize, ...
                binSparseFileHeader.NbrSpots, ...
                mask);

fprintf('matRad: done!\n');

try
    % wait 0.1s for closing all waitbars
    allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar');
    delete(allWaitBarFigures);
    pause(0.1);
catch
end

%% clear all data
delete([MCsquareConfig.CT_File(1:end-4) '.mhd'])
delete([MCsquareConfig.CT_File(1:end-4) '.raw'])
eval(['rmdir ' MCsquareConfig.Output_Directory ' s'])

% cd back
cd(currFolder);

end