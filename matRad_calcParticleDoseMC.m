function dij = matRad_calcParticleDoseMC(ct,stf,pln,cst,nCasePerBixel,calcDoseDirect)
% matRad MCsqaure monte carlo photon dose calculation wrapper
%
% call
%   dij = matRad_calcParticleDoseMc(ct,stf,pln,cst,calcDoseDirect)
%
% input
%   ct:          	matRad ct struct
%   stf:         	atRad steering information struct
%   pln:            matRad plan meta information struct
%   cst:            matRad cst struct
%   nCasePerBixel:  number of histories per beamlet
%   calcDoseDirect: binary switch to enable forward dose
%                   calcualtion
% output
%   dij:            matRad dij struct
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


matRad_cfg = MatRad_Config.instance();

% check if valid machine
if ~strcmp(pln.radiationMode,'protons') || ~strcmp(pln.machine,'generic_MCsquare')
    matRad_cfg.dispError('Wrong radiation modality and/or machine. For now MCsquare requires machine generic_MCsquare!');    
end


if nargin < 5
    % set number of particles simulated per pencil beam
    nCasePerBixel = matRad_cfg.propMC.MCsquare_defaultHistories;
    matRad_cfg.dispInfo('Using default number of Histories per Bixel: %d\n',nCasePerBixel);
end

if nargin < 6
    calcDoseDirect = false;
end

if isfield(pln,'propMC') && isfield(pln.propMC,'outputVariance')
    matRad_cfg.dispWarning('Variance scoring for MCsquare not yet supported.');
end

if ~strcmp(pln.radiationMode,'protons')
    errordlg('MCsquare is only supported for protons');
end

env = matRad_getEnvironment();

%% check if binaries are available
%Executables for simulation
if ispc
    if exist('MCSquare_windows.exe','file') ~= 2
        matRad_cfg.dispError('Could not find MCsquare binary.\n');
    else
        mcSquareBinary = 'MCSquare_windows.exe';
    end
elseif ismac
    if exist('MCsquare_mac','file') ~= 2
        matRad_cfg.dispError('Could not find MCsquare binary.\n');
    else
        mcSquareBinary = './MCsquare_mac';
    end
    %error('MCsquare binaries not available for mac OS.\n');
elseif isunix
    if exist('MCsquare_linux','file') ~= 2
        matRad_cfg.dispError('Could not find MCsquare binary.\n');
    else
        mcSquareBinary = './MCsquare_linux';
    end
end

%Mex interface for import of sparse matrix
if ~calcDoseDirect && ~matRad_checkMexFileExists('matRad_sparseBeamletsReaderMCsquare')
    matRad_cfg.dispWarning('Compiled sparse reader interface not found. Trying to compile it on the fly!');      
    try
        matRad_compileMCsquareSparseReader();        
    catch MException
        matRad_cfg.dispError('Could not find/generate mex interface for reading the sparse matrix. \nCause of error:\n%s\n Please compile it yourself.',MException.message);
    end
end


% set and change to MCsquare binary folder
currFolder = pwd;
fullfilename = mfilename('fullpath');
MCsquareFolder = [fullfilename(1:find(fullfilename==filesep,1,'last')) 'MCsquare' filesep 'bin'];

% cd to MCsquare folder (necessary for binary)
cd(MCsquareFolder);

%Check Materials
if ~exist([MCsquareFolder filesep 'Materials'],'dir') || ~exist(fullfile(MCsquareFolder,'Materials','list.dat'),'file')
    matRad_cfg.dispInfo('First call of MCsquare: unzipping Materials...');    
    unzip('Materials.zip');
    matRad_cfg.dispInfo('Done');
end

% Since MCsquare 1.1 only allows similar resolution in x&y, we do some
% extra checks on that before calling calcDoseInit. First, we make sure a
% dose grid resolution is set in the pln struct
if ~isfield(pln,'propDoseCalc') ...
        || ~isfield(pln.propDoseCalc,'doseGrid') ...
        || ~isfield(pln.propDoseCalc.doseGrid,'resolution') ...
        || ~all(isfield(pln.propDoseCalc.doseGrid.resolution,{'x','y','z'}))
    
    %Take default values
    pln.propDoseCalc.doseGrid.resolution = matRad_cfg.propDoseCalc.defaultResolution;
end

% Now we check for different x/y
if pln.propDoseCalc.doseGrid.resolution.x ~= pln.propDoseCalc.doseGrid.resolution.y
    pln.propDoseCalc.doseGrid.resolution.x = mean([pln.propDoseCalc.doseGrid.resolution.x pln.propDoseCalc.doseGrid.resolution.y]);
    pln.propDoseCalc.doseGrid.resolution.y = pln.propDoseCalc.doseGrid.resolution.x;
    matRad_cfg.dispWarning('Anisotropic resolution in axial plane for dose calculation with MCsquare not possible\nUsing average x = y = %g mm\n',pln.propDoseCalc.doseGrid.resolution.x);
end

%Now we can run calcDoseInit as usual
matRad_calcDoseInit;

%Issue a warning when we have more than 1 scenario
if dij.numOfScenarios ~= 1
    matRad_cfg.dispWarning('MCsquare is only implemented for single scenario use at the moment. Will only use the first Scenario for Monte Carlo calculation!');
end

% prefill ordering of MCsquare bixels
dij.MCsquareCalcOrder = NaN*ones(dij.totalNumOfBixels,1);

% We need to adjust the offset used in matRad_calcDoseInit
mcSquareAddIsoCenterOffset = [dij.doseGrid.resolution.x/2 dij.doseGrid.resolution.y/2 dij.doseGrid.resolution.z/2] ...
                - [dij.ctGrid.resolution.x   dij.ctGrid.resolution.y   dij.ctGrid.resolution.z];
mcSquareAddIsoCenterOffset = mcSquareAddIsoCenterOffset - offset;

% for MCsquare we explicitly downsample the ct to the dose grid (might not
% be necessary in future MCsquare versions with separated grids)
for s = 1:dij.numOfScenarios
    HUcube{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cubeHU{s}, ...
                                dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
end

% Explicitly setting the number of threads for MCsquare, 0 is all available
nbThreads = 0;

% set relative dose cutoff for storage in dose influence matrix, we use the
% default value for the lateral cutoff here
relDoseCutoff = 1 - matRad_cfg.propDoseCalc.defaultLateralCutOff;
% set absolute calibration factor
% convert from eV/g/primary to Gy 1e6 primaries
absCalibrationFactorMC2 = 1.602176e-19 * 1.0e+9;

if isequal(pln.propOpt.bioOptimization,'const_RBExD')
            dij.RBE = 1.1;
            matRad_cfg.dispInfo('matRad: Using a constant RBE of %g\n',dij.RBE);
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



counter = 0;             
for i = 1:length(stf)
    stfMCsquare(i).gantryAngle = mod(180-stf(i).gantryAngle,360); %Different MCsquare geometry
    stfMCsquare(i).couchAngle  = stf(i).couchAngle;
    stfMCsquare(i).isoCenter   = stf(i).isoCenter + mcSquareAddIsoCenterOffset;
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
    matRad_cfg.dispError('Invalid ordering of Beamlets for MCsquare computation!');
end

%% MC computation and dij filling
matRad_writeMCsquareinputAllFiles(MCsquareConfigFile,MCsquareConfig,stfMCsquare);

% run MCsquare
[status,cmdout] = system([mcSquareBinary ' ' MCsquareConfigFile],'-echo');
    
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

matRad_cfg.dispInfo('matRad: done!\n');

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

%For Octave temporarily disable confirmation for recursive rmdir
if strcmp(env,'OCTAVE')    
    rmdirConfirmState = confirm_recursive_rmdir(0);
end
rmdir(MCsquareConfig.Output_Directory,'s');

%Reset to old confirmatoin state
if strcmp(env,'OCTAVE')
    confirm_recursive_rmdir(rmdirConfirmState);
end

% cd back
cd(currFolder);

end
