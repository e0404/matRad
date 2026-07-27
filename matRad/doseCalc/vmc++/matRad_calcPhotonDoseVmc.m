function dij = matRad_calcPhotonDoseVmc(ct, stf, pln, cst, calcDoseDirect, dij)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad vmc++ photon dose calculation wrapper
%
% call
%   dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst,calcDoseDirect,dij)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   calcDoseDirect:             boolian switch to bypass dose influence matrix
%                               computation and directly calculate dose; only makes
%                               sense in combination with matRad_calcDoseDirect.m%
%   dij:                        initialized current-format dij (optional)
% output
%   dij:                        matRad dij struct
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default: dose influence matrix computation
if nargin < 5
    calcDoseDirect = false;
end
if nargin < 6
    dij = [];
end

% Retain the functional entry point while allowing the DoseEngine class to
% supply the current dij schema produced by initDoseCalc.
[dij, ct, numOfVoxels] = matRad_prepareVmcDij(dij, ct, stf);
dij.radiationMode = pln.radiationMode;

% set output level. 0 = no vmc specific output. 1 = print to matlab cmd.
% 2 = open in terminal(s)
verbose = 0;

VMCPath = DoseEngines.matRad_PhotonVmcEngine.getVmcRoot();
if ~isdeployed % only if _not_ running as standalone
    addpath(VMCPath);
end

% check if full dose influence data is required
if calcDoseDirect
    numOfColumnsDij           = length(stf);
    numOfBixelsContainer = 1;
else
    numOfColumnsDij           = dij.totalNumOfBixels;
    numOfBixelsContainer = ceil(dij.totalNumOfBixels / 10);
end

% set up arrays for book keeping
dij.bixelNum = NaN * ones(numOfColumnsDij, 1);
dij.rayNum   = NaN * ones(numOfColumnsDij, 1);
dij.beamNum  = NaN * ones(numOfColumnsDij, 1);

bixelNum = NaN * ones(dij.totalNumOfBixels, 1);
rayNum   = NaN * ones(dij.totalNumOfBixels, 1);
beamNum  = NaN * ones(dij.totalNumOfBixels, 1);

doseTmpContainer        = cell(numOfBixelsContainer, dij.numOfScenarios);
doseTmpContainerError   = cell(numOfBixelsContainer, dij.numOfScenarios);

% Allocate space for dij.physicalDose sparse matrix
for i = 1:dij.numOfScenarios
    dij.physicalDose{i} = spalloc(numOfVoxels, numOfColumnsDij, 1);
    dij.physicalDoseError{i} = spalloc(numOfVoxels, numOfColumnsDij, 1);
end

% set environment variables for vmc++
if ~isfolder(fullfile(VMCPath, 'bin'))
    error('matRad:vmc:EnvironmentNotFound', ...
          'Could not locate the VMC++ environment at %s.', VMCPath);
else
    switch pln.propDoseCalc.vmcOptions.version
        case 'Carleton'
            runsPath    = fullfile(VMCPath, 'run');
        case 'dkfz'
            runsPath    = fullfile(VMCPath, 'runs');
    end
    phantomPath = fullfile(runsPath, 'phantoms');

    setenv('vmc_home', VMCPath);
    setenv('vmc_dir', runsPath);
    setenv('xvmc_dir', VMCPath);

    if isunix
        system(['chmod a+x ' VMCPath filesep 'bin' filesep 'vmc_Linux.exe']);
    end

end

% set consistent random seed (enables reproducibility)
rng(0);

% get default vmc options
VmcOptions = matRad_vmcOptions(pln, ct);

% export CT cube as binary file for vmc++
matRad_exportCtVmc(ct, fullfile(phantomPath, 'matRad_CT.ct'));

% take only voxels inside patient
V = [cst{:, 4}];
V = unique(vertcat(V{:}));

writeCounter        = 0;
readCounter         = 0;
maxNumOfParMCSim    = 0;

% initialize waitbar
matRad_cfg = MatRad_Config.instance();
figureWait = [];
if ~matRad_cfg.disableGUI
    figureWait = waitbar(0, 'calculate dose influence matrix for photons (vmc++)...');
    % show busy state
    set(figureWait, 'pointer', 'watch');
end

fprintf('matRad: VMC++ photon dose calculation...\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams % loop over all beams

    fprintf('Beam %d of %d ...', i, dij.numOfBeams);

    % remember beam and bixel number
    if calcDoseDirect
        dij.beamNum(i)    = i;
        dij.rayNum(i)     = i;
        dij.bixelNum(i)   = i;
    end

    if strcmp(pln.propDoseCalc.vmcOptions.source, 'phsp')
        % set angle-specific vmc++ parameters

        % phsp starts off pointed in the +z direction, with source at -z
        % phsp source gets translated, then rotated (-z, +y, -x) around
        % 0, then pushed to isocenter

        % correct for the source to collimator distance and change units mm -> cm
        translationWorld = stf(i).isoCenter + ...
                           [0 0 pln.propDoseCalc.vmcOptions.SCD + stf(i).sourcePoint_bev(2)];
        translation = DoseEngines.matRad_PhotonVmcEngine.worldToVmcCoordinates(translationWorld, ct) / 10;

        % enter in isocentre
        isocenter = DoseEngines.matRad_PhotonVmcEngine.worldToVmcCoordinates(stf(i).isoCenter, ct) / 10;

        % determine vmc++ rotation angles from gantry and couch
        % angles
        angles = matRad_matRad2vmcSourceAngles(stf(i).gantryAngle, stf(i).couchAngle);

        % set vmc++ parameters
        VmcOptions.source.translation   = translation;
        VmcOptions.source.isocenter     = isocenter;
        VmcOptions.source.angles        = angles;
    end

    % use beam-specific CT name
    VmcOptions.geometry.XyzGeometry.CtFile  = strrep(fullfile(runsPath, 'phantoms', 'matRad_CT.ct'), '\', '/'); % path of density matrix (only needed if input method is 'CT-PHANTOM')

    for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!

        writeCounter = writeCounter + 1;

        % create different seeds for every bixel
        VmcOptions.McControl.rngSeeds = [randi(30000), randi(30000)];

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
        switch pln.propDoseCalc.vmcOptions.source
            case 'beamlet'
                % a) change coordinate system (Isocenter cs-> physical cs) and units mm -> cm
                rayCornersWorld = stf(i).ray(j).rayCorners_SCD + repmat(stf(i).isoCenter, 4, 1);
                rayCorners = DoseEngines.matRad_PhotonVmcEngine.worldToVmcCoordinates(rayCornersWorld, ct) / 10;
                rayCorner1 = rayCorners(1, :);
                rayCorner2 = rayCorners(2, :);
                rayCorner3 = rayCorners(3, :); % vmc needs only three corners (counter-clockwise)
                beamSourceWorld = stf(i).sourcePoint + stf(i).isoCenter;
                beamSource = DoseEngines.matRad_PhotonVmcEngine.worldToVmcCoordinates(beamSourceWorld, ct) / 10;

                % b) swap x and y (CT-standard = [y,x,z])
                rayCorner1 = rayCorner1([2, 1, 3]);
                rayCorner2 = rayCorner2([2, 1, 3]);
                rayCorner3 = rayCorner3([2, 1, 3]);
                beamSource = beamSource([2, 1, 3]);

                % c) set vmc++ parameters
                VmcOptions.source.monoEnergy                    = stf(i).ray(j).energy;                 % photon energy
                % VmcOptions.source.monoEnergy                   = []                  ;                  % use photon spectrum
                VmcOptions.source.beamletEdges                  = [rayCorner1, rayCorner2, rayCorner3];    % counter-clockwise beamlet edges
                VmcOptions.source.virtualPointSourcePosition    = beamSource;                            % virtual beam source position

            case 'phsp'
                % use ray-specific file name for the phsp source (bixelized
                % phsp)
                VmcOptions.source.file_name     = strrep(stf(i).ray(j).phspFileName, '\', '/');
        end

        %% create input file with vmc++ parameters
        outfile = ['MCpencilbeam_temp_', num2str(mod(writeCounter - 1, VmcOptions.run.numOfParMCSim) + 1)];
        matRad_createVmcInput(VmcOptions, fullfile(runsPath, [outfile, '.vmc']));

        % parallelization: only run this block for every numOfParallelMCSimulations!!!
        if mod(writeCounter, VmcOptions.run.numOfParMCSim) == 0 || writeCounter == dij.totalNumOfBixels

            % create batch file (enables parallel processes)
            if writeCounter == dij.totalNumOfBixels && mod(writeCounter, VmcOptions.run.numOfParMCSim) ~= 0
                currNumOfParMCSim = mod(writeCounter, VmcOptions.run.numOfParMCSim);
            else
                currNumOfParMCSim = VmcOptions.run.numOfParMCSim;
            end
            matRad_createVmcBatchFile(currNumOfParMCSim, fullfile(VMCPath, 'run_parallel_simulations.bat'), verbose);

            % save max number of executed parallel simulations
            if currNumOfParMCSim > maxNumOfParMCSim
                maxNumOfParMCSim = currNumOfParMCSim;
            end

            %% perform vmc++ simulation
            currentFolder = pwd;
            folderCleanup = onCleanup(@() cd(currentFolder));
            cd(VMCPath);
            if verbose > 0 % only show output if verbose level > 0
                dos('.\run_parallel_simulations.bat');
                fprintf(['Completed ' num2str(writeCounter) ' of ' num2str(dij.totalNumOfBixels) ' beamlets...\n']);
            else
                [~, ~] = dos('.\run_parallel_simulations.bat');
            end
            clear folderCleanup;

            for k = 1:currNumOfParMCSim
                readCounter     = readCounter + 1;

                % update waitbar
                if ishandle(figureWait)
                    waitbar(writeCounter / dij.totalNumOfBixels, figureWait);
                end

                %% import calculated dose
                idx = regexp(outfile, '_');
                switch pln.propDoseCalc.vmcOptions.version
                    case 'Carleton'
                        filename = sprintf('%s%d.dos', outfile(1:idx(2)), k);
                    case 'dkfz'
                        filename = sprintf('%s%d_%s.dos', outfile(1:idx(2)), k, VmcOptions.scoringOptions.outputOptions.name);
                end
                [bixelDose, bixelDoseError] = matRad_readDoseVmc(fullfile(runsPath, filename), VmcOptions);

                %{
                %%% Don't do any sampling, since the correct error is
                difficult to figure out. We also don't really need it on
                the Graham cluster.

                if ~calcDoseDirect
                    % if not calculating dose directly, sample dose

                    % determine cutoff
                    doseCutoff          = VmcOptions.run.relDoseCutoff*max(bixelDose);

                    % determine which voxels to sample
                    indSample = bixelDose < doseCutoff & bixelDose ~= 0;
                    r = rand(nnz(indSample),1);

                    % sample them
                    thresRand = bixelDose(indSample)./doseCutoff;
                    indKeepSampled = r < thresRand;
                    indKeep = indSample;
                    indKeep(indKeep) = indKeepSampled;

                    bixelDose(indKeep) = doseCutoff;
                    bixelDose(indSample & ~indKeep) = 0;

                end
                %}

                % apply absolute calibration factor
                bixelDoseError  = sqrt((VmcOptions.run.absCalibrationFactorVmc .* bixelDoseError).^2 + (bixelDose .* VmcOptions.run.absCalibrationFactorVmc_err).^2);
                bixelDose       = bixelDose * VmcOptions.run.absCalibrationFactorVmc;

                % Save dose for every bixel in cell array
                doseTmpContainer{mod(readCounter - 1, numOfBixelsContainer) + 1, 1}       = sparse(V, 1, bixelDose(V), numOfVoxels, 1);
                doseTmpContainerError{mod(readCounter - 1, numOfBixelsContainer) + 1, 1}  = sparse(V, 1, bixelDoseError(V), numOfVoxels, 1);

                % save computation time and memory by sequentially filling the
                % sparse matrix dose.dij from the cell array
                if mod(readCounter, numOfBixelsContainer) == 0 || readCounter == dij.totalNumOfBixels
                    if calcDoseDirect
                        if isfield(stf(beamNum(readCounter)).ray(rayNum(readCounter)), 'weight')
                            % score physical dose
                            rayWeight = stf(beamNum(readCounter)).ray(rayNum(readCounter)).weight;
                            if iscell(rayWeight)
                                rayWeight = rayWeight{1};
                            else
                                rayWeight = rayWeight(1);
                            end
                            dij.physicalDose{1}(:, i) = dij.physicalDose{1}(:, i) + ...
                                                        rayWeight * doseTmpContainer{1, 1};
                            dij.physicalDoseError{1}(:, i) = sqrt(dij.physicalDoseError{1}(:, i).^2 + ...
                                                                  (rayWeight * doseTmpContainerError{1, 1}).^2);
                        else
                            error(['No weight available for beam ' num2str(beamNum(readCounter)) ', ray ' num2str(rayNum(readCounter))]);
                        end
                    else
                        % fill entire dose influence matrix
                        dij.physicalDose{1}(:, (ceil(readCounter / numOfBixelsContainer) - 1) * numOfBixelsContainer + 1:readCounter) = ...
                            [doseTmpContainer{1:mod(readCounter - 1, numOfBixelsContainer) + 1, 1}];

                        dij.physicalDoseError{1}(:, (ceil(readCounter / numOfBixelsContainer) - 1) * numOfBixelsContainer + 1:readCounter) = ...
                            [doseTmpContainerError{1:mod(readCounter - 1, numOfBixelsContainer) + 1, 1}];
                    end
                end
            end

        end

    end

    fprintf('Done!\n');
end

%% delete temporary files
delete(fullfile(VMCPath, 'run_parallel_simulations.bat')); % batch file
delete(fullfile(phantomPath, 'matRad_CT.ct'));             % phantom file
for j = 1:maxNumOfParMCSim
    delete(fullfile(runsPath, ['MCpencilbeam_temp_', num2str(mod(j - 1, VmcOptions.run.numOfParMCSim) + 1), '.vmc'])); % vmc inputfile
    switch pln.propDoseCalc.vmcOptions.version
        case 'Carleton'
            filename = sprintf('%s%d.dos', 'MCpencilbeam_temp_', j);
        case 'dkfz'
            filename = sprintf('%s%d_%s.dos', 'MCpencilbeam_temp_', j, VmcOptions.scoringOptions.outputOptions.name);
    end
    delete(fullfile(runsPath, filename));    % vmc outputfile
end

if ishandle(figureWait)
    delete(figureWait);
end

end

function [dij, ct, numOfVoxels] = matRad_prepareVmcDij(dij, ct, stf)
% Initialize legacy functional calls and validate the VMC++ grid contract.

ct = matRad_getWorldAxes(ct);
if isempty(dij)
    dij.ctGrid.resolution = ct.resolution;
    dij.ctGrid.x = ct.x;
    dij.ctGrid.y = ct.y;
    dij.ctGrid.z = ct.z;
    dij.ctGrid.dimensions = ct.cubeDim;
    dij.ctGrid.numOfVoxels = prod(ct.cubeDim);
    dij.doseGrid = dij.ctGrid;
    dij.doseGrid.cubeCoordOffset = [0 0 0];
    dij.numOfBeams = numel(stf);
    dij.numOfScenarios = 1;
    dij.numOfRaysPerBeam = [stf(:).numOfRays];
    dij.totalNumOfBixels = sum([stf(:).totalNumOfBixels]);
    dij.totalNumOfRays = sum(dij.numOfRaysPerBeam);
    dij.weightToMU = 100;
    dij.scaleFactor = 1;
end

numOfVoxels = prod(ct.cubeDim);
if dij.doseGrid.numOfVoxels ~= numOfVoxels || ...
        ~isequal(dij.doseGrid.dimensions, ct.cubeDim)
    error('matRad:vmc:InvalidDoseGrid', 'VMC++ requires the dose grid to match the CT grid.');
end
if dij.numOfScenarios ~= 1 || ct.numOfCtScen ~= 1
    error('matRad:vmc:UnsupportedScenarios', 'VMC++ currently supports only a single nominal CT scenario.');
end

end
