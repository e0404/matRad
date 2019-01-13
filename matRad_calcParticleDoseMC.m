function dij = matRad_calcParticleDoseMC(ct,stf,pln,cst,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad MCsqure monte carlo photon dose calculation wrapper
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

% disable visualiazation by default
if nargin < 5
    visBool = false;
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
VdoseGrid = find(interp3(dij.ctGrid.y,  dij.ctGrid.x,   dij.ctGrid.z,tmpCube, ...
                         dij.doseGrid.y,dij.doseGrid.x',dij.doseGrid.z,'nearest'));

% ignore densities outside of contours
eraseCtDensMask = ones(dij.ctGrid.numOfVoxels,1);
eraseCtDensMask(VctGrid) = 0;
for i = 1:ct.numOfCtScen
    ct.cubeHU{i}(eraseCtDensMask == 1) = -1024;
end

% downsample ct
for s = 1:dij.numOfScenarios
    HUcube{s} =  interp3(dij.ctGrid.y,  dij.ctGrid.x',  dij.ctGrid.z,ct.cubeHU{s}, ...
                         dij.doseGrid.y,dij.doseGrid.x',dij.doseGrid.z,'linear');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbThreads = 4;
nbParallel = 8;%nproc();
if nbParallel > 4
  nbParallel = 2*ceil(nbParallel/nbThreads);
end

% set consistent random seed (enables reproducibility)
%rng(0);
rand('state',0)

% set number of particles simulated per pencil beam
nCasePerBixel = 1000;

% set relative dose cutoff for storage in dose influence matrix
relDoseCutoff = 10^(-3);
% set absolute calibration factor
% convert from sum of dose in Gy for all histories to Gy per 1e6 primaries
absCalibrationFactorMC2 = 1.602176e-19 * nCasePerBixel * 1.0e+6;

% Allocate space for dij.physicalDose sparse matrix
for i = 1:dij.numOfScenarios
    dij.physicalDose{i} = spalloc(dij.doseGrid.numOfVoxels,dij.totalNumOfBixels,1);
end

% helper function for energy selection
round2 = @(a,b)round(a*10^b)/10^b;

% Allocate memory for dose_temp cell array
numOfBixelsContainer = nbParallel;

if ~strcmp(pln.radiationMode,'protons')
    errordlg('MCsquare is only supported for protons');
end

doseTmpContainer = cell(numOfBixelsContainer,dij.numOfScenarios);
if isequal(pln.propOpt.bioOptimization,'const_RBExD')
            dij.RBE = 1.1;
            fprintf(['matRad: Using a constant RBE of 1.1 \n']);
end

% MCsquare settings
MCsquareConfigFile = 'MCsquareConfig.txt';

MCsquareConfig = MatRad_MCsquareConfig;

MCsquareConfig.BDL_Plan_File = 'currBixel.txt';
MCsquareConfig.CT_File       = 'MC2patientCT.mhd';
MCsquareConfig.Num_Threads   = nbThreads;
MCsquareConfig.RNG_Seed      = 1234;
MCsquareConfig.Num_Primaries = nCasePerBixel;

% write patient data
MCsquareBinCubeResolution = [dij.doseGrid.resolution.y ...
                             dij.doseGrid.resolution.x ...
                             dij.doseGrid.resolution.z];   
matRad_writeMhd(HUcube{1},MCsquareBinCubeResolution,MCsquareConfig.CT_File);

% prepare steering for MCsquare
counter = 1;
for i = 1:length(stf)
    for j = 1:stf(i).numOfRays
        for k = 1:stf(i).numOfBixelsPerRay(j)
            
            dij.bixelNum(counter) = k;
            dij.rayNum(counter)   = j;
            dij.beamNum(counter)  = i;
            
            MCsquareStf(counter).gantryAngle = mod(180-stf(i).gantryAngle,360);
            MCsquareStf(counter).couchAngle  = stf(i).couchAngle;
            MCsquareStf(counter).isoCenter   = stf(i).isoCenter-MCsquareBinCubeResolution/2;
            MCsquareStf(counter).energy      = stf(i).ray(j).energy(k);
            MCsquareStf(counter).targetPoint = [-stf(i).ray(j).rayPos_bev(1), stf(i).ray(j).rayPos_bev(3)];
            
            counter = counter + 1;
        end
    end
end

% dont ask me why but we need the materials folder in the current folder
copyfile MC2/Materials/ Materials/

%% MC computation and dij filling
for i = 1:numOfBixelsContainer:dij.totalNumOfBixels
  
    nbRuns = numOfBixelsContainer;
    if (i+nbRuns)> dij.totalNumOfBixels
        nbRuns = dij.totalNumOfBixels-i+1;
    end
    
    for runNb = i:i+nbRuns-1
        
        matRad_writeMCsquareinputFiles(MCsquareConfigFile,MCsquareConfig,MCsquareStf(i));

        [status,cmdout] = system(['MCSquare_windows.exe ' MCsquareConfigFile],'-echo');
        
        % Save dose for every bixel in cell array
        bixelDose = matRad_readMhd([pwd filesep MCsquareConfig.Output_Directory],'Dose.mhd');
        
        % apply relative dose cutoff
        doseCutoff                        = relDoseCutoff*max(bixelDose(:));
        bixelDose(bixelDose < doseCutoff) = 0;

        % apply absolute calibration factor
        bixelDose = bixelDose*absCalibrationFactorMC2;

        sparseDose = sparse(VdoseGrid,1,double(bixelDose(VdoseGrid)),dij.doseGrid.numOfVoxels,1);
        
        doseTmpContainer{mod(runNb-1,numOfBixelsContainer)+1,1} = sparseDose;
        
        delete(MCsquareConfig.BDL_Plan_File)
        delete MCsquareConfig.txt
        
    end
    
    
    
    % save computation time and memory by sequentially filling the
    % sparse matrix dose.dij from the cell array
    dij.physicalDose{1}(:,i:i+nbRuns-1) = [doseTmpContainer{1:nbRuns,1}];

    % Display progress and update text only 200 times
    if mod(i+nbRuns-1,max(1,round(dij.totalNumOfBixels/200))) == 0
        matRad_progress((i+nbRuns-1)/max(1,round(dij.totalNumOfBixels/200)),...
                        floor(dij.totalNumOfBixels/max(1,round(dij.totalNumOfBixels/200))));
        fprintf('\n');          
    end
    
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
delete MC2patientCT.mhd
delete MC2patientCT.raw
rmdir MCsquareOutput s
rmdir Materials s

end