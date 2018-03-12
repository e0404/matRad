function dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst,nCasePerBixel,numOfParallelMCSimulations,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad vmc++ photon dose calculation wrapper
% 
% call
%   dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst,nCasePerBixel,numOfParallelMCSimulations)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   nCasePerBixel:              number of photons simulated per bixel
%   numOfParallelMCSimulations: number of simultaneously performed simulations (optional) 
%   calcDoseDirect:             boolian switch to bypass dose influence matrix
%                               computation and directly calculate dose; only makes
%                               sense in combination with matRad_calcDoseDirect.m%
% output
%   dij:                        matRad dij struct
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default: dose influence matrix computation
if ~exist('calcDoseDirect','var')
    calcDoseDirect = false;
end

% set output level. 0 = no vmc specific output. 1 = print to matlab cmd.
% 2 = open in terminal(s)
verbose = 1;

if ~isdeployed % only if _not_ running as standalone    
    % add path for optimization functions
    matRadRootDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(matRadRootDir,'vmc++'))
end

% meta information for dij
dij.numOfBeams         = pln.propStf.numOfBeams;
dij.numOfVoxels        = prod(ct.cubeDim);
dij.resolution         = ct.resolution;
dij.dimensions         = ct.cubeDim;
dij.numOfScenarios     = 1;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);

% check if full dose influence data is required
if calcDoseDirect 
    numOfColumnsDij           = length(stf);
    numOfBixelsContainer = 1;
else
    numOfColumnsDij           = dij.totalNumOfBixels;
    numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
end

% set up arrays for book keeping
dij.bixelNum = NaN*ones(numOfColumnsDij,1);
dij.rayNum   = NaN*ones(numOfColumnsDij,1);
dij.beamNum  = NaN*ones(numOfColumnsDij,1);

bixelNum = NaN*ones(dij.totalNumOfBixels,1);
rayNum   = NaN*ones(dij.totalNumOfBixels,1);
beamNum  = NaN*ones(dij.totalNumOfBixels,1);

doseTmpContainer = cell(numOfBixelsContainer,dij.numOfScenarios);

% Allocate space for dij.physicalDose sparse matrix
for i = 1:dij.numOfScenarios
    dij.physicalDose{i} = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
end

% set environment variables for vmc++
if exist(['vmc++' filesep 'bin'],'dir') ~= 7
    error(['Could not locate vmc++ environment. ' ...
          'Please provide the files in the correct folder structure at matRadroot' filesep 'vmc++.']);
else
    VMCPath     = fullfile(pwd , 'vmc++');
    runsPath    = fullfile(VMCPath, 'runs');
    phantomPath = fullfile(VMCPath, 'phantoms');

    setenv('vmc_home',VMCPath);
    setenv('vmc_dir',runsPath);
    setenv('xvmc_dir',VMCPath);
    
    if isunix
        system(['chmod a+x ' VMCPath filesep 'bin' filesep 'vmc_Linux.exe']);
    end
    
end

% set consistent random seed (enables reproducibility)
rng(0);

% set number of photons simulated per bixel and number of parallel MC simulations if not specified by user
if nargin < 5
    nCasePerBixel              = 5000;
    if ispc
        numOfParallelMCSimulations = 4;
    elseif isunix
        numOfParallelMCSimulations = 1;
    end
    
    warning(['Number of photons simulated per bixel (nCasePerBixel) and number of parallel MC simulations (numOfParallelMCSimulations) not specified by user. ',...
             'Use default settings with nCasePerBixel = ',num2str(nCasePerBixel),...
             ' and numOfParallelMCSimulations = ',num2str(numOfParallelMCSimulations),...
             ' in vmc++ calculations.'])
elseif nargin < 6
    if ispc
        numOfParallelMCSimulations = 4;
    elseif isunix
        numOfParallelMCSimulations = 1;
    end
    
    warning(['Number of parallel MC simulations (numOfParallelMCSimulations) not specified by user. ',...
             'Use default settings with numOfParallelMCSimulations = ',num2str(numOfParallelMCSimulations),...
             ' in vmc++ calculations.'])    
elseif isunix
    if numOfParallelMCSimulations > 1
        numOfParallelMCSimulations = 1;
    end
    warning(['Running Unix environment: Number of parallel MC simulations (numOfParallelMCSimulations) set to default settings with numOfParallelMCSimulations = ',num2str(numOfParallelMCSimulations),...
             ' in vmc++ calculations.'])    
    
end
    
% set relative dose cutoff for storage in dose influence matrix
relDoseCutoff = 10^(-3);

% set absolute calibration factor
% CALCULATION
% absolute_calibration_factor = 1/D(depth = 100,5mm) -> D(depth = 100,5mm) = 1Gy
% SETUP
% SAD = 1000mm, SCD = 500mm, bixelWidth = 5mm, IC = [240mm,240mm,240mm]
% fieldsize@IC = 105mm x 105mm, phantomsize = 81 x 81 x 81 = 243mm x 243mm x 243mm
% rel_Dose_cutoff = 10^(-3), ncase = 500000/bixel
absCalibrationFactorVmc = 99.818252282632300;

% set general vmc++ parameters
% 1 source
VmcOptions.beamletSource.myName       = 'source 1';                        % name of source
VmcOptions.beamletSource.monitorUnits = 1;                                 
VmcOptions.beamletSource.spectrum     = ['./spectra/var_6MV.spectrum'];    % energy spectrum source (only used if no mono-Energy given)
VmcOptions.beamletSource.charge       = 0;                                 % charge (-1,0,1)
% 2 transport parameter
VmcOptions.McParameter.automatic_parameter = 'yes';                        % if yes, automatic transport parameters are used
% 3 MC control
VmcOptions.McControl.ncase  = nCasePerBixel;                               % number of histories
VmcOptions.McControl.nbatch = 10;                                          % number of batches
% 4 variance reduction
VmcOptions.varianceReduction.repeatHistory      = 0.251;
VmcOptions.varianceReduction.splitPhotons       = 'yes';   
VmcOptions.varianceReduction.photonSplitFactor = -40;  
% 5 quasi random numbers
VmcOptions.quasi.base      = 2;                                                 
VmcOptions.quasi.dimension = 60;                                             
VmcOptions.quasi.skip      = 1;                                              
% 6 geometry
VmcOptions.geometry.XyzGeometry.methodOfInput = 'CT-PHANTOM';              % input method ('CT-PHANTOM', 'individual', 'groups') 
VmcOptions.geometry.XyzGeometry.Ct            = 'CT';                      % name of geometry
VmcOptions.geometry.XyzGeometry.CtFile        = ['./phantoms/matRad_CT.ct']; % path of density matrix (only needed if input method is 'CT-PHANTOM')
% 7 scoring manager
VmcOptions.scoringOptions.startInGeometry               = 'CT';            % geometry in which partciles start their transport
VmcOptions.scoringOptions.doseOptions.scoreInGeometries = 'CT';            % geometry in which dose is recorded
VmcOptions.scoringOptions.doseOptions.scoreDoseToWater  = 'yes';           % if yes output is dose to water
VmcOptions.scoringOptions.outputOptions.name            = 'CT';            % geometry for which dose output is created (geometry has to be scored)
VmcOptions.scoringOptions.outputOptions.dumpDose        = 2;               % output format (1: format=float, Dose + deltaDose; 2: format=short int, Dose)

% export CT cube as binary file for vmc++
matRad_exportCtVmc(ct, fullfile(phantomPath, 'matRad_CT.ct'));

% take only voxels inside patient
V = [cst{:,4}];
V = unique(vertcat(V{:}));

writeCounter                  = 0;
readCounter                   = 0;
maxNumOfParallelMcSimulations = 0;

% initialize waitbar
figureWait = waitbar(0,'VMC++ photon dose influence matrix calculation..');

fprintf('matRad: VMC++ photon dose calculation... ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams % loop over all beams
       
   % remember beam and bixel number
    if calcDoseDirect
        dij.beamNum(i)    = i;
        dij.rayNum(i)     = i;
        dij.bixelNum(i)   = i;
    end
    
    for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
        
        writeCounter = writeCounter + 1;

        % create different seeds for every bixel
        VmcOptions.McControl.rngSeeds = [randi(30000),randi(30000)];

        % remember beam and bixel number
        if ~calcDoseDirect
           dij.beamNum(writeCounter)  = i;
           dij.rayNum(writeCounter)   = j;
           dij.bixelNum(writeCounter) = j;
        end
        beamNum(writeCounter)  = i;
        rayNum(writeCounter)   = j;
        bixelNum(writeCounter) = j;
        
        % set ray specific vmc++ parameters
        % a) change coordinate system (Isocenter cs-> physical cs) and units mm -> cm
        rayCorner1 = (stf(i).ray(j).rayCorners_SCD(1,:) + stf(i).isoCenter)/10;              
        rayCorner2 = (stf(i).ray(j).rayCorners_SCD(2,:) + stf(i).isoCenter)/10;
        rayCorner3 = (stf(i).ray(j).rayCorners_SCD(3,:) + stf(i).isoCenter)/10; %vmc needs only three corners (counter-clockwise)
        beamSource = (stf(i).sourcePoint + stf(i).isoCenter)/10;
        
        % b) swap x and y (CT-standard = [y,x,z])
        rayCorner1 = rayCorner1([2,1,3]);              
        rayCorner2 = rayCorner2([2,1,3]);
        rayCorner3 = rayCorner3([2,1,3]);
        beamSource  = beamSource([2,1,3]);
        
        % c) set vmc++ parameters
        VmcOptions.beamletSource.monoEnergy                = stf(i).ray(j).energy;                 % photon energy
        %VmcOptions.beamletSource.monoEnergy                 = []                  ;                  % use photon spectrum
        VmcOptions.beamletSource.beamletEdges               = [rayCorner1,rayCorner2,rayCorner3];    % counter-clockwise beamlet edges
        VmcOptions.beamletSource.virtualPointSourcePosition = beamSource;                            % virtual beam source position
        
        % create inputfile with vmc++ parameters
        outfile = ['MCpencilbeam_temp_',num2str(mod(writeCounter-1,numOfParallelMCSimulations)+1)];
        matRad_createVmcInput(VmcOptions,fullfile(runsPath, [outfile,'.vmc']));
        
        % parallelization: only run this block for every numOfParallelMCSimulations!!!
        if mod(writeCounter,numOfParallelMCSimulations) == 0 || writeCounter == dij.totalNumOfBixels
            
            % create batch file (enables parallel processes)
            if writeCounter == dij.totalNumOfBixels && mod(writeCounter,numOfParallelMCSimulations) ~= 0
                currNumOfParallelMcSimulations = mod(writeCounter,numOfParallelMCSimulations);
            else
                currNumOfParallelMcSimulations = numOfParallelMCSimulations;
            end
            matRad_createVmcBatchFile(currNumOfParallelMcSimulations,fullfile(VMCPath,'run_parallel_simulations.bat'),verbose);
            
            % save max number of executed parallel simulations
            if currNumOfParallelMcSimulations > maxNumOfParallelMcSimulations 
                maxNumOfParallelMcSimulations = currNumOfParallelMcSimulations;
            end
            
            % perform vmc++ simulation
            current = pwd;
            cd(VMCPath);
            if verbose > 0 % only show output if verbose level > 0
                dos('run_parallel_simulations.bat');
                fprintf(['Completed ' num2str(writeCounter) ' of ' num2str(dij.totalNumOfBixels) ' beamlets...\n']);
            else
                [dummyOut1,dummyOut2] = dos('run_parallel_simulations.bat'); % supress output by assigning dummy output arguments
            end
            cd(current);
            
            for k = 1:currNumOfParallelMcSimulations
                readCounter = readCounter + 1;
                
                % Display progress
                if verbose == 0
                   % matRad_progress(readCounter,dij.totalNumOfBixels);
                end
                
                % update waitbar
                waitbar(writeCounter/dij.totalNumOfBixels);
                
                % import calculated dose
                idx = regexp(outfile,'_');
                [bixelDose,~] = matRad_readDoseVmc(fullfile(VMCPath, 'runs',...
                                                     [outfile(1:idx(2)),num2str(k), '_', VmcOptions.scoringOptions.outputOptions.name, '.dos']));

                % apply relative dose cutoff
                doseCutoff                        = relDoseCutoff*max(bixelDose);
                bixelDose(bixelDose < doseCutoff) = 0;

                % apply absolute calibration factor
                bixelDose = bixelDose*absCalibrationFactorVmc;

                % Save dose for every bixel in cell array
                doseTmpContainer{mod(readCounter-1,numOfBixelsContainer)+1,1} = sparse(V,1,bixelDose(V),dij.numOfVoxels,1);
                
                % save computation time and memory by sequentially filling the 
                % sparse matrix dose.dij from the cell array
                if mod(readCounter,numOfBixelsContainer) == 0 || readCounter == dij.totalNumOfBixels
                    if calcDoseDirect
                        if isfield(stf(beamNum(readCounter)).ray(rayNum(readCounter)),'weight')
                            % score physical dose
                            dij.physicalDose{1}(:,i) = dij.physicalDose{1}(:,i) + stf(beamNum(readCounter)).ray(rayNum(readCounter)).weight * doseTmpContainer{1,1};
                        else
                            error(['No weight available for beam ' num2str(beamNum(readCounter)) ', ray ' num2str(rayNum(readCounter))]);
                        end
                    else
                        % fill entire dose influence matrix
                        dij.physicalDose{1}(:,(ceil(readCounter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:readCounter) = ...
                            [doseTmpContainer{1:mod(readCounter-1,numOfBixelsContainer)+1,1}];
                    end
                end
            end
            
        end
        
    end
end

% delete temporary files
delete(fullfile(VMCPath, 'run_parallel_simulations.bat')); % batch file
delete(fullfile(phantomPath, 'matRad_CT.ct'));             % phantom file
for j = 1:maxNumOfParallelMcSimulations
    delete(fullfile(runsPath, ['MCpencilbeam_temp_',num2str(mod(j-1,numOfParallelMCSimulations)+1),'.vmc'])); % vmc inputfile
    delete(fullfile(runsPath, ['MCpencilbeam_temp_',num2str(mod(j-1,numOfParallelMCSimulations)+1),'_',...
                                    VmcOptions.scoringOptions.doseOptions.scoreInGeometries,'.dos']));    % vmc outputfile
end

try
  % wait 0.1s for closing all waitbars
  allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar'); 
  delete(allWaitBarFigures);
  pause(0.1); 
catch
end
