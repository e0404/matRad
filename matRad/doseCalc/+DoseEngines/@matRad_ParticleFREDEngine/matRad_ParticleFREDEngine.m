classdef matRad_ParticleFREDEngine < DoseEngines.matRad_MonteCarloEngineAbstract
    % Engine for particle dose calculation using FRED Monte Carlo algorithm
    % for more information see superclass
    % DoseEngines.matRad_MonteCarloEngineAbstract
    %
    %
    % The following parameters for the FRED engine can be tuned by the user. In
    % order to do so, specify the desired value in: pln.propDoseCalc.
    % [s]: string/character array
    % [b]: boolean
    % [i]: integer
    % [f]: float/double/any non strictly integer number
    %
    %
    % HUclamping:            [b] allows for clamping of HU table. Default: true
    % HUtable:               [s] HU table name. Example: 'internal', 'matRad_default_FRED'
    % externalCalculation    [b/s] off (default): run FRED
    %                              t/'write'  : Only write simulation parameter files
    %                              'path'     : read simulation files from 'path'
    %
    % sourceModel            [s] see AvailableSourceModels, {'gaussian', 'emittance', 'sigmaSqrModel'}
    % useGPU                 [b] trigger use of GPU (if available)
    % roomMaterial           [s] material of the patient surroundings. Example:
    %                            'vacuum', 'Air'
    % printOutput            [b] 't: FRED output is mirrored to Matlab console, f: no output is printed'
    % numHistoriesDirect     [i]
    % numHistoriesPerBeamlet [i]
    % scorers                [c] cell array with specified scorers. Example:
    %                            'Dose', 'LETd'
    % primaryMass            [f] mass of the primary ion (in Da). Default value for
    %                             protons: 1.0727
    % numOfNucleons          [i] number of nucleons. Default for protons: 1

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2023 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Constant)
        possibleRadiationModes = {'protons'}
        name                   = 'FRED'
        shortName              = 'FRED'
    end

    properties (SetAccess = protected, GetAccess = public)

        defaultHUtable        = 'matRad_default_FredMaterialConverter'
        availableSourceModels = {'gaussian', 'emittance', 'sigmaSqrModel'}
        calcBioDose
        currentVersion
        availableVersions = {'3.70.0'}  % Or higher.
        radiationMode

    end

    properties
        externalCalculation = 'write'
        useGPU = true
        calcLET = false
        constantRBE
        HUclamping
        scorers
        HUtable
        sourceModel
        roomMaterial
        printOutput
        primaryMass
        numOfNucleons
        ignoreOutsideDensities
        workingDir
        forceDijFormatVersion
    end

    properties (Dependent)
        dijFormatVersion
    end

    properties (SetAccess = private, Hidden)
        patientFilename      = 'CTpatient.mhd'
        runInputFilename     = 'fred.inp'
        regionsFilename      = 'regions.inp'
        funcsFilename        = 'funcs.inp'
        planFilename         = 'plan.inp'
        fieldsFilename       = 'fields.inp'
        layersFilename       = 'layers.inp'
        beamletsFilename     = 'beamlets.inp'
        planDeliveryFilename = 'planDelivery.inp'

        hLutLimits = [-1000, 1375]   % Default FRED values

        conversionFactor = 1e6      % Used to scale the FRED dose to matRad normalization
        MCrunFolder
        inputFolder
        regionsFolder
        planFolder

        HUcube
    end

    methods

        function this = matRad_ParticleFREDEngine(pln)
            % Constructor
            %
            % call
            %   engine = DoseEngines.matRad_DoseEngineFRED(ct,stf,pln,cst)
            %

            matRad_cfg = MatRad_Config.instance();
            if nargin < 1
                pln = [];
            end

            % call superclass constructor
            this = this@DoseEngines.matRad_MonteCarloEngineAbstract(pln);

            if nargin > 0
                if isfield(pln, 'radiationMode')
                    this.radiationMode = pln.radiationMode;
                end
            end

            if isempty(this.workingDir)
                this.workingDir = fullfile(matRad_cfg.primaryUserFolder, 'FRED');
            end

            if ~exist(this.workingDir, 'dir')
                mkdir(this.workingDir);
                matRad_cfg.dispWarning('FRED root folder not found, this should not happen!');
            end
        end

    end

    methods (Access = protected)

        dij = calcDose(this, ct, cst, stf)

        function dij = initDoseCalc(this, ct, cst, stf)

            matRad_cfg = MatRad_Config.instance();

            dij = initDoseCalc@DoseEngines.matRad_MonteCarloEngineAbstract(this, ct, cst, stf);

            dij = this.allocateQuantityMatrixContainers(dij, {'physicalDose'});

            % Issue a warning when we have more than 1 scenario
            if dij.numOfScenarios ~= 1
                matRad_cfg.dispWarning( ...
                                       ['FRED is only implemented for single scenario use at the moment. '...
                                        'Will only use the first Scenario for Monte Carlo calculation!'] ...
                                      );
            end

            % Check for model consistency
            if ~isempty(this.bioModel) && isa(this.bioModel, 'matRad_LQLETbasedModel')
                this.calcBioDose = true;
            else
                this.calcBioDose = 0;
            end

            % Limit RBE calculation to proton models for the time being
            if this.calcBioDose

                switch this.radiationMode

                    case 'protons'

                        dij = this.loadBiologicalData(cst, dij);
                        dij = this.allocateQuantityMatrixContainers(dij, {'mAlphaDose', 'mSqrtBetaDose'});

                        % Only considering LET based models
                        this.calcLET = true;
                    otherwise
                        matRad_cfg.dispWarning('biological dose calculation not supported for radiation modality: %s', this.radiationMode);
                        this.calcBioDose = false;
                end
            end

            if isa(this.bioModel, 'matRad_ConstantRBE')
                dij.RBE = this.bioModel.RBE;
            end

            if this.calcLET
                this.scorers = [this.scorers, {'LETd'}];
                % Allocate containers for both LET*Dose and dose weighted
                % LET. This last is used for biological calculation as well
                dij = this.allocateQuantityMatrixContainers(dij, {'mLETDose', 'mLETd'});
            end

        end

        function writeTreeDirectory(this)

            % Loop over the scenarios
            if this.multScen.totNumScen > 1
                for scenIdx = 1:this.multScen.totNumScen
                    this.writeTreeDirectoryForRun(scenIdx);
                end
            else

                this.writeTreeDirectoryForRun(0);
            end
        end

        function writeTreeDirectoryForRun(this, scenIdx)

            [~, runFolderName] = fileparts(this.MCrunFolder);
            if scenIdx == 0
                tailRun = '';
            else
                tailRun = sprintf('_%d', scenIdx);
            end

            folderName = sprintf('%s%s', this.MCrunFolder, tailRun);
            if ~exist(folderName, 'dir')
                mkdir(folderName);
            end

            % write input folder
            folderName = strrep(this.inputFolder, runFolderName, sprintf('%s%s', runFolderName, tailRun));
            if ~exist(folderName, 'dir')
                mkdir(folderName);
            end

            % build MCrun/inp/regions
            folderName = strrep(this.regionsFolder, runFolderName, sprintf('%s%s', runFolderName, tailRun));
            if ~exist(folderName, 'dir')
                mkdir(folderName);
            end

            % build MCrun/inp/plan
            folderName = strrep(this.planFolder, runFolderName, sprintf('%s%s', runFolderName, tailRun));
            if ~exist(folderName, 'dir')
                mkdir(folderName);
            end
        end

        %% Write files functions

        writeRunFile(~, fName)

        writeRegionsFile(this, fName, stf)

        writePlanDeliveryFile(this, fName, stf)

        writePlanFile(this, fName, stf, scenIdx)

        function writeFredInputAllFiles(this, stf)

            % write fred.inp file
            for scenIdx = 1:this.multScen.totNumScen
                [~, runFolderName] = fileparts(this.MCrunFolder);

                if this.multScen.totNumScen > 1
                    tailRun = sprintf('_%d', scenIdx);
                else
                    tailRun = '';
                end

                runFilename = fullfile(strrep(this.MCrunFolder, runFolderName, sprintf('%s%s', runFolderName, tailRun)), this.runInputFilename);
                this.writeRunFile(runFilename);

                % write region/region.inp file
                regionFilename = fullfile(strrep(this.regionsFolder, runFolderName, sprintf('%s%s', runFolderName, tailRun)), this.regionsFilename);
                this.writeRegionsFile(regionFilename);

                if ~strcmp(this.HUtable, 'internal')
                    hlutFilename = fullfile(strrep(this.regionsFolder, runFolderName, sprintf('%s%s', runFolderName, tailRun)), 'hLut.inp');
                    this.writeHlutFile(hlutFilename, scenIdx);
                end

                % write plan file
                planFile = fullfile(strrep(this.planFolder, runFolderName, sprintf('%s%s', runFolderName, tailRun)), this.planFilename);
                this.writePlanFile(planFile, stf, scenIdx);

                % write planDelivery file
                planDeliveryFile = fullfile(strrep(this.planFolder, runFolderName, sprintf('%s%s', runFolderName, tailRun)), ...
                                            this.planDeliveryFilename);
                this.writePlanDeliveryFile(planDeliveryFile);
            end
        end

        function writeCTs(this)

            patientMetadata.imageOrigin = [0 0 0];
            patientMetadata.resolution  = [this.doseGrid.resolution.x, this.doseGrid.resolution.y, this.doseGrid.resolution.z];
            patientMetadata.datatype = 'int16';

            if this.multScen.totNumScen > 1
                [~, runFolderName] = fileparts(this.MCrunFolder);

                for scenIdx = 1:this.multScen.totNumScen
                    fileNamePatient = fullfile(strrep(this.regionsFolder, runFolderName, sprintf('%s_%d', runFolderName, scenIdx)), ...
                                               this.patientFilename);
                    ctIdx = this.multScen.linearMask(scenIdx, 1);
                    matRad_writeMHD(fileNamePatient, this.HUcube{ctIdx}, patientMetadata);
                end

            else
                fileNamePatient = fullfile(this.regionsFolder, this.patientFilename);
                matRad_writeMHD(fileNamePatient, this.HUcube{1}, patientMetadata);
            end

        end

        function writeHlutFile(this, fileName, scenIdx)

            matRad_cfg = MatRad_Config.instance();

            if ~exist('scenIdx', 'var') || isempty(scenIdx)
                scenIdx = 1;
            end

            mainFolder        = fullfile(matRad_cfg.matRadSrcRoot, 'hluts');
            userDefinedFolder = fullfile(matRad_cfg.primaryUserFolder, 'hluts');
            fredDefinedFolder = fullfile(matRad_cfg.matRadSrcRoot, 'doseCalc', 'FRED', 'hluts');

            % Collect all the subfolders

            searchPath = [strsplit(genpath(mainFolder), pathsep)'; ...
                          strsplit(genpath(userDefinedFolder), pathsep)'; ...
                          strsplit(genpath(fredDefinedFolder), pathsep)'];

            searchPath(cellfun(@isempty, searchPath)) = [];

            % Check for existence of folder paths
            searchPath = searchPath(cellfun(@isfolder, searchPath));

            availableHLUTs = cellfun(@(x) dir([x, filesep, '*.txt']), searchPath, 'UniformOutput', false);
            availableHLUTs = cell2mat(availableHLUTs);

            hLUTindex = find(strcmp([this.HUtable, '.txt'], {availableHLUTs.name}));

            if isempty(hLUTindex)
                errString = sprintf('Cannot open hLut: %s. Available hLut files are: ', hLutFile);
                errString = [errString, sprintf('\ninternal')];
                for hLUTindex = 1:numel(availableHLUTs)
                    errString = [errString, sprintf('\n%s', strrep(fullfile(availableHLUTs(hLUTindex).folder, ...
                                                                            availableHLUTs(hLUTindex).name), '\', '\\'))];
                end
                matRad_cfg.dispError(errString);

            else
                selectedHlutfile = fullfile(availableHLUTs(hLUTindex).folder, availableHLUTs(hLUTindex).name);
                selectedHlut = this.readHlutFileToStruct(selectedHlutfile);
            end

            % Apply relative range shift
            selectedHlut.RSP = selectedHlut.RSP * (1 + this.multScen.relRangeShift(scenIdx));

            % How do we handle absolute shift? Do we insert a slab of water
            % in front of the patient?

            % Write the hlut back
            this.writeHlutFileFromStruct(fileName, selectedHlut);

        end

        function materials = readHlutFileToStruct(this, fileName)

            matRad_cfg = MatRad_Config.instance();
            fid = fopen(fileName, 'r');
            if fid == -1
                matRad_cfg.dispError('Cannot open file.');
            end

            % --- Read header line ---
            headerLine = fgetl(fid);

            % Remove "matColumns:" and split column names
            headerLine = strrep(headerLine, 'matColumns:', '');
            colNames = strsplit(strtrim(headerLine));

            % --- Read data lines ---
            data = [];
            while ~feof(fid)
                line = strtrim(fgetl(fid));
                if startsWith(line, 'mat:')
                    line = strrep(line, 'mat:', '');
                    values = sscanf(line, '%f')';
                    data = [data; values];
                end
            end

            fclose(fid);

            % --- Assign main fields ---
            materials = struct();

            materials.HU   = data(:, strcmp(colNames, 'HU'))';
            materials.rho  = data(:, strcmp(colNames, 'rho'))';
            materials.RSP  = data(:, strcmp(colNames, 'RSP'))';
            materials.Ipot = data(:, strcmp(colNames, 'Ipot'))';
            materials.Lrad = data(:, strcmp(colNames, 'Lrad'))';

            % --- Composition fields ---
            compStartIdx = find(strcmp(colNames, 'C'));  % first composition column
            compNames = colNames(compStartIdx:end);

            materials.materialComposition = struct();

            for i = 1:length(compNames)
                materials.materialComposition.(compNames{i}) = ...
                    data(:, compStartIdx + i - 1)';
            end

        end

        function writeHlutFileFromStruct(this, fileName, materials)

            matRad_cfg = MatRad_Config.instance();

            fid = fopen(fileName, 'w');
            if fid == -1
                matRad_cfg.dispError('Cannot open file.');
            end

            % --- Main fields ---
            mainFields = {'HU', 'rho', 'RSP', 'Ipot', 'Lrad'};

            % --- Composition fields ---
            compFields = fieldnames(materials.materialComposition);

            % --- Write header ---
            fprintf(fid, 'matColumns: ');
            fprintf(fid, '%s ', mainFields{:});
            fprintf(fid, '%s ', compFields{:});
            fprintf(fid, '\n');

            % --- Number of materials ---
            n = numel(materials.HU);

            % --- Write rows ---
            for i = 1:n
                fprintf(fid, 'mat: ');

                % Write main properties
                for k = 1:numel(mainFields)
                    value = materials.(mainFields{k})(i);
                    fprintf(fid, '%g ', value);
                end

                % Write composition values
                for k = 1:numel(compFields)
                    value = materials.materialComposition.(compFields{k})(i);
                    fprintf(fid, '%g ', value);
                end

                fprintf(fid, '\n');
            end

            fclose(fid);

        end

        function dij = loadBiologicalData(this, cst, dij)
            matRad_cfg = MatRad_Config.instance();

            matRad_cfg.dispInfo('Initializing biological dose calculation...\n');

            dij.ax              = zeros(dij.doseGrid.numOfVoxels, 1);
            dij.bx              = zeros(dij.doseGrid.numOfVoxels, 1);

            cstDownsampled = matRad_setOverlapPriorities(cst);

            % resizing cst to dose cube resolution
            cstDownsampled = matRad_resizeCstToGrid(cstDownsampled, dij.ctGrid.x, dij.ctGrid.y, dij.ctGrid.z, ...
                                                    dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);

            % retrieve photon LQM parameter for the current dose grid voxels
            [dij.ax, dij.bx] = matRad_getPhotonLQMParameters(cstDownsampled, dij.doseGrid.numOfVoxels, this.VdoseGrid);

        end

        function dij = allocateQuantityMatrixContainers(this, dij, names)

            % if this.calcDoseDirect
            %     numOfBixelsContainer = 1;
            % else
            %     numOfBixelsContainer = dij.totalNumOfBixels;
            % end

            % Loop over all requested quantities
            for n = 1:numel(names)

                dij.(names{n}) = cell(size(this.multScen.scenMask));

                % Now preallocate a matrix in each active scenario using the
                % scenmask
                if this.calcDoseDirect
                    dij.(names{n})(this.multScen.scenMask) = {zeros(dij.doseGrid.numOfVoxels, this.numOfColumnsDij)};
                else
                    % We preallocate a sparse matrix with sparsity of
                    % 1e-3 to make the filling slightly faster
                    % TODO: the preallocation could probably
                    % have more accurate estimates

                    dij.(names{n})(this.multScen.scenMask) = {spalloc(dij.doseGrid.numOfVoxels, ...
                                                                      this.numOfColumnsDij, ...
                                                                      round(prod(dij.doseGrid.numOfVoxels, ...
                                                                                 this.numOfColumnsDij) * 1e-3)) ...
                                                             };
                end
            end

        end

    end

    methods

        function setDefaults(this)
            setDefaults@DoseEngines.matRad_MonteCarloEngineAbstract(this);

            matRad_cfg = MatRad_Config.instance();

            this.HUclamping             = false;
            this.HUtable                = this.defaultHUtable;
            this.externalCalculation    = 'off';
            this.useGPU                 = false;
            this.sourceModel            = this.availableSourceModels{1};
            this.roomMaterial           = 'Air';
            this.printOutput            = true;
            this.numHistoriesDirect     = matRad_cfg.defaults.propDoseCalc.numHistoriesDirect;
            this.numHistoriesPerBeamlet = matRad_cfg.defaults.propDoseCalc.numHistoriesPerBeamlet;
            this.scorers                = {'Dose'};
            this.primaryMass            = 1.00727;
            this.numOfNucleons          = 1;
            this.outputMCvariance       = false;
            this.constantRBE            = NaN;
            this.ignoreOutsideDensities = false;

        end

        function writeHlut(this, hLutFile, fileName)

            matRad_cfg = MatRad_Config.instance();

            mainFolder        = fullfile(matRad_cfg.matRadSrcRoot, 'hluts');
            userDefinedFolder = fullfile(matRad_cfg.primaryUserFolder, 'hluts');
            fredDefinedFolder = fullfile(matRad_cfg.matRadSrcRoot, 'doseCalc', 'FRED', 'hluts');

            % Collect all the subfolders

            searchPath = [strsplit(genpath(mainFolder), pathsep)'; ...
                          strsplit(genpath(userDefinedFolder), pathsep)'; ...
                          strsplit(genpath(fredDefinedFolder), pathsep)'];

            searchPath(cellfun(@isempty, searchPath)) = [];

            % Check for existence of folder paths
            searchPath = searchPath(cellfun(@isfolder, searchPath));

            availableHLUTs = cellfun(@(x) dir([x, filesep, '*.txt']), searchPath, 'UniformOutput', false);
            availableHLUTs = cell2mat(availableHLUTs);

            hLUTindex = find(strcmp([hLutFile, '.txt'], {availableHLUTs.name}));

            if ~isempty(hLUTindex)
                selectedHlutfile = fullfile(availableHLUTs(hLUTindex).folder, availableHLUTs(hLUTindex).name);

                template = fileread(selectedHlutfile);

                newLut = fopen(fileName, 'w');
                fprintf(newLut, template);
                fclose(newLut);
            else

                errString = sprintf('Cannot open hLut: %s. Available hLut files are: ', hLutFile);
                errString = [errString, sprintf('\ninternal')];
                for hLUTindex = 1:numel(availableHLUTs)

                    errString = [errString, sprintf('\n%s', ...
                                                    strrep(fullfile(availableHLUTs(hLUTindex).folder, availableHLUTs(hLUTindex).name), '\', '\\'))];

                end
                matRad_cfg.dispError(errString);
            end

        end

    end

    methods (Static)

        function cmdString = cmdCall(newCmdString)
            persistent fredCmdCall
            if nargin > 0
                fredCmdCall = newCmdString;
            elseif isempty(fredCmdCall)
                if ispc
                    %                    fredCmdCall = 'wsl if [ -f ~/.fredenv.sh ] ; then source ~/.fredenv.sh ; fi; fred';
                    fredCmdCall = 'fred ';
                elseif isunix
                    fredCmdCall = 'if [ -f ~/.fredenv.sh ] ; then source ~/.fredenv.sh ; fi; fred';
                else
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('OS not supported for FRED!');
                end
            end
            cmdString = fredCmdCall;
        end

        function availableVersions = getAvailableVersions()
            % Function to get available FRED version

            matRad_cfg = MatRad_Config.instance();

            currCmdCall = DoseEngines.matRad_ParticleFREDEngine.cmdCall;

            availableVersions = [];

            [status, cmdOut] = system([currCmdCall, ' -listVers']);

            if status == 0
                nLidx = regexp(cmdOut, '\n') + 6; % 6 because of tab
                nVersions = numel(nLidx) - 1;

                for versIdx = 1:nVersions
                    availableVersions = [availableVersions, {cmdOut(nLidx(versIdx):nLidx(versIdx) + 5)}];
                end

            else
                matRad_cfg.dispError('Something wrong occurred in checking FRED available version. Please check correct FRED installation');
            end
        end

        function [available, msg] = isAvailable(pln, machine)
            % see superclass for information

            msg = [];
            available = false;

            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            % checkBasic
            try
                checkBasic = isfield(machine, 'meta') && isfield(machine, 'data');

                % check modality
                checkModality = any(strcmp(DoseEngines.matRad_ParticleFREDEngine.possibleRadiationModes, machine.meta.radiationMode));

                preCheck = checkBasic && checkModality;

                if ~preCheck
                    return
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return
            end
            available = preCheck;
        end

        function execCheck = checkExec()

            matRad_cfg = MatRad_Config.instance();

            % Check if I can obtain FRED version
            try
                ver = DoseEngines.matRad_ParticleFREDEngine.getVersion();
                if ~isempty(ver)
                    execCheck = true;
                else
                    execCheck = false;

                    msg = sprintf(['Couldn''t call FRED executable. ' ...
                                   'Please set the correct path with DoseEngines.matRad_ParticleFREDEngine.cmdCall(''path/to/executable''). ' ...
                                   'Current value is ''%s'''], DoseEngines.matRad_ParticleFREDEngine.cmdCall);

                    matRad_cfg.dispError(msg);
                end
            catch
                execCheck = false;

                msg = sprintf(['Couldn''t call FRED executable. ' ...
                               'Please set the correct path with DoseEngines.matRad_ParticleFREDEngine.cmdCall(''path/to/executable''). ' ...
                               'Current value is ''%s'''], DoseEngines.matRad_ParticleFREDEngine.cmdCall);
                matRad_cfg.dispError(msg);
            end
        end

        function version = getVersion()
            % Function to get current default FRED version
            matRad_cfg = MatRad_Config.instance();

            try
                [status, cmdOut] = system([DoseEngines.matRad_ParticleFREDEngine.cmdCall, ' -vn']);

                if status == 0
                    % Extract the version number using a regular expression
                    versionMatch = regexp(cmdOut, '\d+\.\d+\.\d+', 'match', 'once');
                    if ~isempty(versionMatch)
                        version = versionMatch;
                    else
                        matRad_cfg.dispError('Unable to parse FRED version from output. Please check FRED installation.');
                        version = [];
                    end
                else
                    matRad_cfg.dispError('Error occurred while checking FRED installation. Please check FRED installation.');
                    version = [];
                end
            catch
                matRad_cfg.dispWarning('Something wrong occurred in checking FRED installation. Please check correct FRED installation');
                version = [];
            end
        end

        % end

        function dijMats = readSparseDataV2(fID, fileFormatVersion, cubeDim, numberOfBixels, nComponents)
            matRad_cfg = MatRad_Config.instance();

            values = zeros([nComponents, 0]);
            voxelIndices = [];
            colIndices = [];

            for i = 1:numberOfBixels
                % Read Beamlet
                bixNum = fread(fID, 1, "int32");
                numVox  = fread(fID, 1, "int32");

                colIndices(end + 1:end + numVox) = i;
                currVoxelIndices = fread(fID, numVox, "uint32") + 1;
                values(:, end + 1:end + numVox) = fread(fID, [nComponents numVox], "float32");

                % x and y components have been permuted in CT
                [indY, indX, indZ] = ind2sub(cubeDim, currVoxelIndices);

                voxelIndices(end + 1:end + numVox) = sub2ind(cubeDim([2, 1, 3]), indX, indY, indZ);
                if matRad_cfg.logLevel > 2
                    matRad_progress(i, numberOfBixels);
                end
            end
            for c = 1:nComponents
                dijMats{c} = sparse(voxelIndices, colIndices, values(c, :), prod(cubeDim), numberOfBixels);
            end
        end

        function dijMats = readSparseDataV3(fID, fileFormatVersion, cubeDim, numberOfBixels, nComponents)
            matRad_cfg = MatRad_Config.instance();

            allBixelMeta = fread(fID, [3 numberOfBixels], "uint32"); % (#bixels, PBidx, FID, PBID)

            % Size information for each component
            componentDataSize = fread(f, nComponents, "uint32");

            % Data
            dijMatrices  = cell(1, nComponents);
            for c = 1:nComponents
                PBidxs       = fread(fID, componentDataSize(c), "uint32") + 1;
                voxelIndices = fread(fID, componentDataSize(c), "uint32") + 1;
                values       = fread(fID, componentDataSize(c), "float32");

                % x and y components have been permuted in CT
                [indY, indX, indZ] = ind2sub(cubeDim, voxelIndices);
                voxelIndices = sub2ind(cubeDim([2, 1, 3]), indX, indY, indZ);

                dijMatrices{c} = sparse(voxelIndices, PBidxs, values, prod(cubeDim), numberOfBixels);
            end
        end

        function dijMatrices = readSparseDijBin(fName)
            % FRED function to read sparseDij in .bin format
            % call
            %   readSparseDijBin(fName)
            %
            % input
            %   fName: filename to read
            %
            % output
            %   dijMatrix: dij structure
            matRad_cfg = MatRad_Config.instance();

            f = fopen(fName, 'r', 'l');

            try
                % Header
                fileFormatVersion = fread(f, 1, "int32");
                dims = fread(f, 3, "int32");
                res = fread(f, 3, "float32");
                offset = fread(f, 3, "float32");
                if fileFormatVersion > 20
                    orientation = fread(f, 9, "float32");
                end
                nComponents = fread(f, 1, "int32");
                numberOfBixels = fread(f, 1, "int32");

                matRad_cfg.dispInfo('Reading FRED dij with %d components of size %dx%d (voxels x beamlets) with cubeDim = %dx%dx%d\n', ...
                                    nComponents, prod(dims), numberOfBixels, dims(1), dims(2), dims(3));

                if fileFormatVersion < 30
                    dijMatrices = DoseEngines.matRad_ParticleFREDEngine.readSparseDataV2(f, fileFormatVersion, dims, numberOfBixels, nComponents);
                else
                    dijMatrices = DoseEngines.matRad_ParticleFREDEngine.readSparseDataV3(f, fileFormatVersion, dims, numberOfBixels, nComponents);
                end

                fclose(f);
            catch ME
                fclose(f);
                matRad_cfg.dispError('unable to load file %s: %s', fName, ME.message);
            end
        end

        % Used to check against a machine file if a specific quantity can be
        % computed.
        function q = providedQuantities(machine)
            q = {'physicalDose', 'LET'};
        end

        [doseCubeV, letdCubeV, fileName] = readSimulationOutput(runFolder, calcDoseDirect, calcLET)

    end

    methods (Access = private)

        function updatePaths(obj, rootFolder)

            if ~strcmp(rootFolder, obj.workingDir)
                obj.workingDir  = rootFolder;
            end

            obj.MCrunFolder     = fullfile(obj.workingDir, 'MCrun');
            obj.inputFolder     = fullfile(obj.MCrunFolder, 'inp');
            obj.regionsFolder   = fullfile(obj.inputFolder, 'regions');
            obj.planFolder      = fullfile(obj.inputFolder, 'plan');

        end

        function [radiationMode] = updateRadiationMode(this, value)

            % This function also resets the values for primary mass and number
            % of nucleons. Used for possible future extension to multiple
            % ion species
            matRad_cfg = MatRad_Config.instance();

            if any(strcmp(value, this.possibleRadiationModes))
                radiationMode = value;
            else
                matRad_cfg.dispError('Invalid radiation modality: %s', value);
            end

            switch radiationMode
                case 'protons'
                    this.primaryMass = 1.00727; % Da
                    this.numOfNucleons = 1;
                    matRad_cfg.dispInfo('Default values for priamry mass and number of nucleons set.');
                otherwise
                    matRad_cfg.dispError('Only proton dose calculation available with this version of FRED');

            end

            matRad_cfg.dispWarning('Selected radiation modality: %s with primary mass: %2.3f', radiationMode, this.primaryMass);
        end

        function isLower = isVersionLower(this, version)
            % This function directly looks at FRED installation, not at
            % the current FRED version stored in the class property.
            fredVersion = this.getVersion();

            isLower = false;

            if ~isempty(fredVersion)
                % Decompose the current version for comparison
                vdiff = sscanf(fredVersion, '%d.%d.%d') - sscanf(version, '%d.%d.%d');
                firstdiff = find(vdiff, 1, 'first');
                isLower = ~isempty(firstdiff) && vdiff(firstdiff) < 0;
            end
        end

    end

    methods

        function set.sourceModel(this, value)
            matRad_cfg = MatRad_Config.instance();

            valid = ischar(value) && any(strcmp(value, this.availableSourceModels));

            if valid
                this.sourceModel = value;
            else

                matRad_cfg.dispWarning('Unable to set source model:%s, setting default:%s', value, this.availableSourceModels{1});
                this.sourceModel = this.availableSourceModels{1};

            end

        end

        function version = get.currentVersion(this)

            if isempty(this.currentVersion)
                version = this.getVersion();
                this.currentVersion = version;
            else
                version = this.currentVersion;
            end

        end

        function v = get.dijFormatVersion(this)

            matRad_cfg = MatRad_Config.instance();

            if ~isempty(this.forceDijFormatVersion)
                v = this.forceDijFormatVersion;
            elseif this.isVersionLower('3.76.0')
                v = 20;
            else
                v = 31;
            end

            % FRED version <= 3.70.0 does not allow dij version
            % selection and only works with ifFormatVersion < 21
            if this.isVersionLower('3.76.0')
                if v > 20
                    matRad_cfg.dispWarning('FRED version %s does not support ijFormatVersions>20. Version 20 will be used!');
                    v = 20;
                end
            else
                if v == 20
                    matRad_cfg.dispWarning('FRED version %s does no longer support ijFormatVersions 20. Version 21 will be used!');
                    v = 21;
                end
            end

        end

        function set.radiationMode(this, value)
            if ischar(value)
                if ~isempty(this.radiationMode) && ~strcmp(this.radiationMode, value)
                    this.radiationMode = value;
                    this.updateRadiationMode(this.radiationMode);
                elseif isempty(this.radiationMode)
                    this.radiationMode = value;
                end
            end

        end

        function set.workingDir(obj, pathValue)
            obj.workingDir = pathValue;
            obj.updatePaths(pathValue);
        end

        function set.externalCalculation(this, value)
            % Set exportCalculation value, available options are:
            %  - false:    (default) runs the FRED simulation (requires FRED installation)
            %  - write/1:  triggers the file export
            %  - 'path':   simulation data will be loaded from the specified
            %              path. Full simulation directory path should be provided.
            %              Example: 'matRadRoot/userdata/FRED/'

            if isnumeric(value) || islogical(value)
                switch value
                    case 1
                        this.externalCalculation = 'write';
                    case 0
                        this.externalCalculation = 'off';
                end
            elseif ischar(value)

                if any(strcmp(value, {'write', 'off'}))
                    this.externalCalculation = value;
                elseif isfolder(value)
                    this.externalCalculation = value;

                    this.updatePaths(value);
                end
            end
        end

    end
end
