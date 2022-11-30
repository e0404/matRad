function dij = matRad_calcParticleDoseMCsquare(ct,stf,pln,cst,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad MCsquare Monte Carlo proton dose calculation wrapper
%
% call
%   dij = matRad_calcParticleDoseMCsquare(ct,stf,pln,cst,calcDoseDirect)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct            
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

% Instance of MatRad_Config class
matRad_cfg = MatRad_Config.instance();

% Handle inputs
if nargin < 5
    calcDoseDirect = false;
end

% initialize waitbar
figureWait = waitbar(0,'calculate dose influence matrix with MCsquare...');
% prevent closure of waitbar and show busy state
set(figureWait,'pointer','watch');

% check if valid machine
if ~strcmp(pln.radiationMode,'protons')
    matRad_cfg.dispError('Wrong radiation modality. MCsquare only supports protons!');    
end

% Load class variables in pln
% for calcDoseDirect, this is already done in superior function
if ~isfield(pln,'propMC') || ~isa(pln.propMC,'MatRad_MCsquareConfig')
    pln = matRad_cfg.getDefaultClass(pln,'propMC','MatRad_MCsquareConfig');
end

% load default parameters in case they haven't been set yet
pln = matRad_cfg.getDefaultProperties(pln,'propDoseCalc');

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
    matRad_cfg.dispInfo('Done!\n');
end

% Since MCsquare 1.1 only allows similar resolution in x&y, we do some
% extra checks on that before calling calcDoseInit. First, we make sure a
% dose grid resolution is set in the pln struct
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
relDoseCutOff = 1 - matRad_cfg.propDoseCalc.defaultLateralCutOff;
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

% switch for using existing BDL file (e.g. to fit matRad basedata), 
% or generate BDL file from matRad base data using MCsquareBDL
if isfield(pln,'loadExistingBDL') && ~isempty(pln.loadExistingBDL)
    % use existing BDL file
    bdFile = pln.loadExistingBDL;
    
else
    % fit and create BDL file using selected machine file
    bdFile = [machine.meta.machine '.txt'];

    % override base data in case of APM, it is not needed here
    machine.data = MatRad_HeterogeneityConfig.overrideBaseData(machine.data);

    % Calculate MCsquare base data
    % Argument stf is optional, if given, calculation only for energies given in stf
    MCsquareBDL = MatRad_MCsquareBaseData(machine);

    %matRad_createMCsquareBaseDataFile(bdFile,machine,1);
    MCsquareBDL = MCsquareBDL.writeMCsquareData([MCsquareFolder filesep 'BDL' filesep bdFile]);
    MCsquareBDL = MCsquareBDL.saveMatradMachine('savedMatRadMachine');
    
end

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
                              
                pln.propMC.BDL_Machine_Parameter_File = ['BDL/' bdFile];
                pln.propMC.BDL_Plan_File = 'currBixels.txt';
                pln.propMC.CT_File       = 'MC2patientCT.mhd';
                pln.propMC.Num_Threads   = nbThreads;
                pln.propMC.RNG_Seed      = 1234;              

                % turn simulation of individual beamlets
                pln.propMC.Beamlet_Mode = ~calcDoseDirect;
                % turn of writing of full dose cube
                pln.propMC.Dose_MHD_Output = calcDoseDirect;
                % turn on sparse output
                pln.propMC.Dose_Sparse_Output = ~calcDoseDirect;
                % set threshold of sparse matrix generation
                pln.propMC.Dose_Sparse_Threshold = relDoseCutOff;
                
                %Matrices for LET
                if pln.propDoseCalc.calcLET
                    pln.propMC.LET_MHD_Output		 = calcDoseDirect;
                    pln.propMC.LET_Sparse_Output	 = ~calcDoseDirect;
                end
                 
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
                                         round(stf(i).ray(j).weight(stf(i).ray(j).energy == stfMCsquare(i).energies(k))*pln.propMC.numHistories)];
                                     
                                     totalWeights = totalWeights + stf(i).ray(j).weight(stf(i).ray(j).energy == stfMCsquare(i).energies(k));
                                 else
                                     stfMCsquare(i).energyLayer(k).numOfPrimaries = [stfMCsquare(i).energyLayer(k).numOfPrimaries ...
                                         pln.propMC.numHistories];
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

                %% Write config files
                % override HU_Density_Conversion_File and HU_Material_Conversion_File in case of Heterogeneity density sampling
                if isfield(ct,'modulated') && ct.modulated
                    % copy and override with default HU conversion files
                    copyfile(MCsquareConfig.HU_Density_Conversion_File,'Scanners/densitySampling/HU_Density_Conversion.txt')
                    copyfile(MCsquareConfig.HU_Material_Conversion_File,'Scanners/densitySampling/HU_Material_Conversion.txt')

                    % prepare sampled densities and combine with HU
                    sampledDensities(1,:) = 6000:5999+length(ct.sampledDensities);
                    sampledDensities(2,:) = ct.sampledDensities;

                    % write sampled densities
                    fID = fopen('Scanners/densitySampling/HU_Density_Conversion.txt','a');
                    fprintf(fID,'\n%i       %.3f',sampledDensities);
                    fclose(fID);
                    
                    % write material conversion
                    fID = fopen('Scanners/densitySampling/HU_Material_Conversion.txt','a');
                    fprintf(fID,'\n6000    40      # Schneider_Lung');
%                     fprintf(fID,'\n6000    17      # Water');
                    fclose(fID);

                    % set custom HU conversion files to be used by MCsquare
                    MCsquareConfig.HU_Density_Conversion_File = 'Scanners/densitySampling/HU_Density_Conversion.txt';
                    MCsquareConfig.HU_Material_Conversion_File = 'Scanners/densitySampling/HU_Material_Conversion.txt';
                end

                % write patient data
                MCsquareBinCubeResolution = [dij.doseGrid.resolution.x ...
                    dij.doseGrid.resolution.y ...
                    dij.doseGrid.resolution.z];
                
                pln.propMC.writeMhd(HUcube{ctScen},MCsquareBinCubeResolution);

                % write config file
                pln.propMC.writeMCsquareinputAllFiles(MCsquareConfigFile,stfMCsquare);
                
                %% MC computation and dij filling
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
                        [pln.propMC.Output_Directory filesep 'Sparse_Dose.bin'], ...
                        dij.doseGrid.dimensions, ...
                        dij.totalNumOfBixels, ...
                        mask);
                    
                    %Read sparse LET
                    if pln.propDoseCalc.calcLET
                        dij.mLETDose{1} =  absCalibrationFactorMC2 * matRad_sparseBeamletsReaderMCsquare ( ...
                            [pln.propMC.Output_Directory filesep 'Sparse_LET.bin'], ...
                            dij.doseGrid.dimensions, ...
                            dij.totalNumOfBixels, ...
                            mask);
                    end
                else
                    %Read dose cube
                    cube = pln.propMC.readMhd('Dose.mhd');
                    dij.physicalDose{1} = absCalibrationFactorMC2 * totalWeights * ...
                        sparse(VdoseGrid,ones(numel(VdoseGrid),1), ...
                            cube(VdoseGrid), ...
                            dij.doseGrid.numOfVoxels,1);
                    
                    %Read LET cube
                    if pln.propDoseCalc.calcLET
                        cube = pln.propMC.readMhd('LET.mhd');
                        dij.mLETDose{1} = absCalibrationFactorMC2 * totalWeights * ...
                            sparse(VdoseGrid,ones(numel(VdoseGrid),1), ...
                                cube(VdoseGrid), ...
                                dij.doseGrid.numOfVoxels,1);
                    end

                    % Postprocessing for dij:
                    % This is already the combined dose over all bixels, so all parameters are 1 in this case
                    dij = rmfield(dij,'MCsquareCalcOrder');

                    dij.numOfBeams = 1;
                    dij.beamNum = 1;
                    dij.bixelNum = 1;
                    dij.rayNum = 1;
                    dij.totalNumOfBixels = 1;
                    dij.totalNumOfRays = 1;
                    dij.numOfRaysPerBeam = 1;
                end
                
                % reorder influence matrix to comply with matRad default ordering
                if pln.propMC.Beamlet_Mode
                    dij.physicalDose{1} = dij.physicalDose{1}(:,MCsquareOrder);
                    if pln.propDoseCalc.calcLET
                        dij.mLETDose{1} = dij.mLETDose{1}(:,MCsquareOrder);
                    end
                end
                
                matRad_cfg.dispInfo('Simulation finished!\n');
                
                %% Clear data
                delete([pln.propMC.CT_File(1:end-4) '.*']);
                delete('currBixels.txt');
                delete('MCsquareConfig.txt');

                %For Octave temporarily disable confirmation for recursive rmdir
                if strcmp(env,'OCTAVE')
                    rmdirConfirmState = confirm_recursive_rmdir(0);
                end
                rmdir(pln.propMC.Output_Directory,'s');

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

% Order fields for easier comparison between different dijs
dij = orderfields(dij);

%% cd back
cd(currFolder);

if ishandle(figureWait)
    delete(figureWait);
end

end
