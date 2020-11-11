function dij = matRad_calcParticleDoseMCsquare(ct,stf,pln,cst,nCasePerBixel,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad MCsquare Monte Carlo proton dose calculation wrapper
%
% call
%   dij = matRad_calcParticleDoseMc(ct,stf,pln,cst,calcDoseDirect)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   nCasePerBixel               number of histories per beamlet (nCasePerBixel > 1),
%                               max stat uncertainity (0 < nCasePerBixel < 1)
%   calcDoseDirect:             binary switch to enable forward dose
%                               calculation
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


matRad_cfg = MatRad_Config.instance();

% initialize waitbar
figureWait = waitbar(0,'calculate dose influence matrix with MCsquare...');
% prevent closure of waitbar and show busy state
set(figureWait,'pointer','watch');

% check if valid machine
if ~strcmp(pln.radiationMode,'protons')
    matRad_cfg.dispError('Wrong radiation modality . MCsquare only supports protons!');    
end

if nargin < 5
    % set number of particles simulated per pencil beam
    nCasePerBixel = matRad_cfg.propMC.MCsquare_defaultHistories;
    matRad_cfg.dispInfo('Using default number of Histories per Bixel: %d\n',nCasePerBixel);
end
% switch between either using max stat uncertainity or total number of
% cases
if (nCasePerBixel < 1)
    maxStatUncertainty = true;
else
    maxStatUncertainty = false;
end

if nargin < 6
    calcDoseDirect = false;
end

if isfield(pln,'propMC') && isfield(pln.propMC,'outputVariance')
    matRad_cfg.dispWarning('Variance scoring for MCsquare not yet supported.');
end

if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'calcLET') 
    pln.propDoseCalc.calcLET = matRad_cfg.propDoseCalc.defaultCalcLET;
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
        mcSquareBinary = 'chmod a+x MCsquare_linux && ./MCsquare_linux';
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


% Explicitly setting the number of threads for MCsquare, 0 is all available
nbThreads = 0;

% set relative dose cutoff for storage in dose influence matrix, we use the
% default value for the lateral cutoff here
relDoseCutoff = 1 - matRad_cfg.propDoseCalc.defaultLateralCutOff;
% set absolute calibration factor
% convert from eV/g/primary to Gy 1e6 primaries
absCalibrationFactorMC2 = 1.602176e-19 * 1.0e+9;

% book keeping - this is necessary since pln is not used in optimization or
% matRad_calcCubes
if strcmp(pln.bioParam.model,'constRBE')
    dij.RBE = pln.bioParam.RBE;
end

scenCount = 0;

% for MCsquare we explicitly downsample the ct to the dose grid (might not
% be necessary in future MCsquare versions with separated grids)
for s = 1:ct.numOfCtScen
    HUcube{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cubeHU{s}, ...
        dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
end

%BDL File
bdFile = [machine.meta.machine '.txt'];

% bdFile = 'BDL_matRad.txt'; %use for baseData fit 
MCsquareBDL = MatRad_MCsquareBaseData(machine,stf);
%matRad_createMCsquareBaseDataFile(bdFile,machine,1);
% MCsquareBDL = MCsquareBDL.saveMatradMachine('test');
MCsquareBDL = MCsquareBDL.writeMCsquareData([MCsquareFolder filesep 'BDL' filesep bdFile]);
%movefile(bdFile,[MCsquareFolder filesep 'BDL/' bdFile]);
% MCsquareBDL = MCsquareBDL.saveMatradMachine('testMachine');

for shiftScen = 1:pln.multScen.totNumShiftScen
    
    % manipulate isocenter
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter + pln.multScen.isoShift(shiftScen,:);
    end
    
    for ctScen = 1:pln.multScen.numOfCtScen
        for rangeShiftScen = 1:pln.multScen.totNumRangeScen
            if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                
                %For direct dose calculation
                totalWeights = 0;
                
                %Count the scenarios
                scenCount = scenCount + 1;
                
                % We need to adjust the offset used in matRad_calcDoseInit
                mcSquareAddIsoCenterOffset = [dij.doseGrid.resolution.x/2 dij.doseGrid.resolution.y/2 dij.doseGrid.resolution.z/2] ...
                    - [dij.ctGrid.resolution.x   dij.ctGrid.resolution.y   dij.ctGrid.resolution.z];
                mcSquareAddIsoCenterOffset = mcSquareAddIsoCenterOffset - offset;
                
                % MCsquare settings
                MCsquareConfigFile = sprintf('MCsquareConfig.txt');
                
                MCsquareConfig = MatRad_MCsquareConfig;
                
                MCsquareConfig.BDL_Machine_Parameter_File = ['BDL/' bdFile];
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
                
                %Matrices for LET
                if pln.propDoseCalc.calcLET
                    MCsquareConfig.LET_MHD_Output		 = calcDoseDirect;
                    MCsquareConfig.LET_Sparse_Output	 = ~calcDoseDirect;
                end
                
                
                % write patient data
                MCsquareBinCubeResolution = [dij.doseGrid.resolution.x ...
                    dij.doseGrid.resolution.y ...
                    dij.doseGrid.resolution.z];
                
                matRad_writeMhd(HUcube{ctScen},MCsquareBinCubeResolution,MCsquareConfig.CT_File);
                
                
                
                counter = 0;
                for i = 1:length(stf)
                    %Let's check if we have a unique or no range shifter, because MCsquare
                    %only allows one range shifter type per field which can be IN or OUT
                    %per spot
                    raShiField = [];
                    for j = 1:stf(i).numOfRays
                        if isfield(stf(i).ray(j),'rangeShifter')
                            raShiField = [raShiField stf(i).ray(j).rangeShifter(:).ID];
                        else
                            raShiField = [raShiField zeros(size(stf(i).ray(j).energies))];
                        end
                    end
                    
                    raShiField = unique(raShiField); %unique range shifter
                    raShiField(raShiField == 0) = []; %no range shifter
                    if numel(raShiField) > 1
                        matRad_cfg.dispError('MCsquare does not support different range shifter IDs per field! Aborting.\n');
                    end
                    
                    if ~isempty(raShiField)
                        stfMCsquare(i).rangeShifterID = raShiField;
                        stfMCsquare(i).rangeShifterType = 'binary';
                    else
                        stfMCsquare(i).rangeShifterID = 0;
                        stfMCsquare(i).rangeShifterType = 'binary';
                    end
                    
                    stfMCsquare(i).gantryAngle = mod(180-stf(i).gantryAngle,360); %Different MCsquare geometry
                    stfMCsquare(i).couchAngle  = stf(i).couchAngle;
                    stfMCsquare(i).isoCenter   = stf(i).isoCenter + mcSquareAddIsoCenterOffset;
                    stfMCsquare(i).energies    = unique([stf(i).ray.energy]);
                    stfMCsquare(i).SAD         = stf(i).SAD;
                    
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
                             
                             %Check if ray has a spot in the current energy layer
                             if any(stf(i).ray(j).energy == stfMCsquare(i).energies(k))
                                 %Set up the ray geometries and add current ray to energy
                                 %layer
                                 energyIx = find(stf(i).ray(j).energy == stfMCsquare(i).energies(k));
                                 stfMCsquare(i).energyLayer(k).rayNum   = [stfMCsquare(i).energyLayer(k).rayNum j];
                                 stfMCsquare(i).energyLayer(k).bixelNum = [stfMCsquare(i).energyLayer(k).bixelNum energyIx];
                                 stfMCsquare(i).energyLayer(k).targetPoints = [stfMCsquare(i).energyLayer(k).targetPoints; ...
                                     -stf(i).ray(j).rayPos_bev(1) stf(i).ray(j).rayPos_bev(3)];
                                 
                                 %Number of primaries depending on beamlet-wise or field-based compuation (direct dose calculation)                    
                                 if calcDoseDirect
                                     stfMCsquare(i).energyLayer(k).numOfPrimaries = [stfMCsquare(i).energyLayer(k).numOfPrimaries ...
                                         round(stf(i).ray(j).weight(stf(i).ray(j).energy == stfMCsquare(i).energies(k))*MCsquareConfig.Num_Primaries)];
                                     
                                     totalWeights = totalWeights + stf(i).ray(j).weight(stf(i).ray(j).energy == stfMCsquare(i).energies(k));
                                 else
                                     stfMCsquare(i).energyLayer(k).numOfPrimaries = [stfMCsquare(i).energyLayer(k).numOfPrimaries ...
                                         MCsquareConfig.Num_Primaries];
                                 end
                                 
                                 %Now add the range shifter
                                 raShis = stf(i).ray(j).rangeShifter(energyIx);
                                                                  
                                 %sanity check range shifters
                                 raShiIDs = unique([raShis.ID]);
                                 %raShiIDs = raShiIDs(raShiIDs ~= 0);
                                 
                                 if ~isscalar(raShiIDs)
                                     matRad_cfg.dispError('MCsquare only supports one range shifter setting (on or off) per energy! Aborting.\n');
                                 end
                                 
                                 stfMCsquare(i).energyLayer(k).rangeShifter = raShis(1);
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
                mcSquareCall = [mcSquareBinary ' ' MCsquareConfigFile];
                matRad_cfg.dispInfo(['Calling Monte Carlo Engine: ' mcSquareCall]);
                [status,cmdout] = system(mcSquareCall,'-echo');

                mask = false(dij.doseGrid.numOfVoxels,1);
                mask(VdoseGrid) = true;
                
                % read output
                if ~calcDoseDirect
                    %Read Sparse Matrix
                    dij.physicalDose{1} = absCalibrationFactorMC2 * matRad_sparseBeamletsReaderMCsquare ( ...
                        [MCsquareConfig.Output_Directory filesep 'Sparse_Dose.bin'], ...
                        dij.doseGrid.dimensions, ...
                        dij.totalNumOfBixels, ...
                        mask);
                    
                    %Read sparse LET
                    if pln.propDoseCalc.calcLET
                        dij.mLETDose{1} =  absCalibrationFactorMC2 * matRad_sparseBeamletsReaderMCsquare ( ...
                            [MCsquareConfig.Output_Directory filesep 'Sparse_LET.bin'], ...
                            dij.doseGrid.dimensions, ...
                            dij.totalNumOfBixels, ...
                            mask);
                        
                        dij.MC_tallies{1} = 'LET';
                    end
                else
                    %Read dose cube
                    cube = matRad_readMhd(MCsquareConfig.Output_Directory,'Dose.mhd');
                    dij.physicalDose{1} = absCalibrationFactorMC2 * totalWeights * ...
                        sparse(VdoseGrid,ones(numel(VdoseGrid),1), ...
                            cube(VdoseGrid), ...
                            dij.doseGrid.numOfVoxels,1);
                    
                    %Read LET cube
                    if pln.propDoseCalc.calcLET
                        cube = matRad_readMhd(MCsquareConfig.Output_Directory,'LET.mhd');
                        dij.mLETDose{1} = absCalibrationFactorMC2 * totalWeights * ...
                            sparse(VdoseGrid,ones(numel(VdoseGrid),1), ...
                                cube(VdoseGrid), ...
                                dij.doseGrid.numOfVoxels,1);
                        
                        dij.MC_tallies{1} = 'LET';
                    end
                end
                
                % reorder influence matrix to comply with matRad default ordering
                if MCsquareConfig.Beamlet_Mode
                    dij.physicalDose{1} = dij.physicalDose{1}(:,MCsquareOrder);
                    if pln.propDoseCalc.calcLET
                        dij.mLETDose{1} = dij.mLETDose{1}(:,MCsquareOrder);
                    end
                end
                
                matRad_cfg.dispInfo('Simulation finished!\n');
                
                %% Clear data
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
                
            end
        end
    end
    
    % manipulate isocenter
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter - pln.multScen.isoShift(shiftScen,:);
    end   
end

%% cd back
cd(currFolder);

if ishandle(figureWait)
    delete(figureWait);
end

end
