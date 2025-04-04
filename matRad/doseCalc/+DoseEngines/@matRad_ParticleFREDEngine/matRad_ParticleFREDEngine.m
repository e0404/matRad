classdef matRad_ParticleFREDEngine < DoseEngines.matRad_MonteCarloEngineAbstract
% Engine for particle dose calculation using FRED Monte Carlo algorithm 
% for more informations see superclass
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
%                              t/'write'  : Only write simulation paramter files
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
        possibleRadiationModes = {'protons'};
        name                   = 'FRED';
        shortName              = 'FRED';
    end

    properties (SetAccess = protected, GetAccess = public)
        
        defaultHUtable        = 'matRad_default_FredMaterialConverter';
        AvailableSourceModels = {'gaussian', 'emittance', 'sigmaSqrModel'};
        defaultDijFormatVersion = '20';
      
        calcBioDose;
        currentVersion;
        availableVersions = {'3.70.0'}; % Or higher.
        radiationMode;

    end

    properties
        dijFormatVersion;
        externalCalculation;
        useGPU;
        calcLET;
        constantRBE;
        HUclamping;
        scorers;
        HUtable;
        sourceModel;
        roomMaterial;
        printOutput;
        primaryMass;
        numOfNucleons;
        ignoreOutsideDensities;
    end


    properties (SetAccess = private, Hidden)
        patientFilename      = 'CTpatient.mhd';
        runInputFilename     = 'fred.inp';
        regionsFilename      = 'regions.inp';
        funcsFilename        = 'funcs.inp';
        planFilename         = 'plan.inp';
        fieldsFilename       = 'fields.inp';
        layersFilename       = 'layers.inp';
        beamletsFilename     = 'beamlets.inp';
        planDeliveryFilename = 'planDelivery.inp';

        hLutLimits = [-1000,1375];  % Default FRED values
        
        conversionFactor = 1e6;     % Used to scale the FRED dose to matRad normalization

        FREDrootFolder;

        MCrunFolder;
        inputFolder;
        regionsFolder;
        planFolder;
        dijReaderHandle;
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

            if isempty(this.FREDrootFolder)
                this.FREDrootFolder = fullfile(matRad_cfg.primaryUserFolder, 'FRED');
            end
            
            if ~exist(this.FREDrootFolder, 'dir')
                mkdir(this.FREDrootFolder);
                matRad_cfg.dispWarning('FRED root folder not found, this should not happen!');
            end

        end

    end

    methods(Access = protected)

        dij = calcDose(this,ct,cst,stf)

        function dij = initDoseCalc(this,ct,cst,stf)

            matRad_cfg = MatRad_Config.instance();

            dij = initDoseCalc@DoseEngines.matRad_MonteCarloEngineAbstract(this,ct,cst,stf); 
            
            dij = this.allocateQuantityMatrixContainers(dij, {'physicalDose'});

            %Issue a warning when we have more than 1 scenario
            if dij.numOfScenarios ~= 1
                matRad_cfg.dispWarning('FRED is only implemented for single scenario use at the moment. Will only use the first Scenario for Monte Carlo calculation!');
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

                        dij = this.loadBiologicalData(cst,dij);
                        dij = this.allocateQuantityMatrixContainers(dij,{'mAlphaDose', 'mSqrtBetaDose'});

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

            if ~exist(this.MCrunFolder, 'dir')
                mkdir(this.MCrunFolder);
            end

            % write input folder
            if ~exist(this.inputFolder, 'dir')
                mkdir(this.inputFolder);
            end

            % build MCrun/inp/regions
            if ~exist(this.regionsFolder, 'dir')
                mkdir(this.regionsFolder);
            end

            % build MCrun/inp/plan
            if ~exist(this.planFolder, 'dir')
                mkdir(this.planFolder);
            end
        end


        %% Write files functions

        writeRunFile(~, fName)
                
        writeRegionsFile(this,fName, stf)

        writePlanDeliveryFile(this, fName, stf)
  
        writePlanFile(this,fName, stf)

        function writeFredInputAllFiles(this,stf)
    
            %write fred.inp file
            runFilename = fullfile(this.MCrunFolder, this.runInputFilename);
            this.writeRunFile(runFilename);
        
            %write region/region.inp file
            regionFilename = fullfile(this.regionsFolder, this.regionsFilename);
            this.writeRegionsFile(regionFilename);
        
            %write plan file
            planFile = fullfile(this.planFolder, this.planFilename);
            this.writePlanFile(planFile,stf);

            %write planDelivery file

            planDeliveryFile = fullfile(this.planFolder,this.planDeliveryFilename);
            this.writePlanDeliveryFile(planDeliveryFile);
        end

        function dij = loadBiologicalData(this,cst,dij)
           matRad_cfg = MatRad_Config.instance();
            
            matRad_cfg.dispInfo('Initializing biological dose calculation...\n');
            
            dij.ax              = zeros(dij.doseGrid.numOfVoxels,1);
            dij.bx              = zeros(dij.doseGrid.numOfVoxels,1);
            
            cstDownsampled = matRad_setOverlapPriorities(cst);
            
            % resizing cst to dose cube resolution
            cstDownsampled = matRad_resizeCstToGrid(cstDownsampled,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,...
                dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);
            
            % retrieve photon LQM parameter for the current dose grid voxels
            [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cstDownsampled,dij.doseGrid.numOfVoxels,this.VdoseGrid);
            
        end


        function dij = allocateQuantityMatrixContainers(this,dij,names)

            % if this.calcDoseDirect
            %     numOfBixelsContainer = 1;
            % else
            %     numOfBixelsContainer = dij.totalNumOfBixels;
            % end
            
            %Loop over all requested quantities
            for n = 1:numel(names)

                dij.(names{n}) = cell(size(this.multScen.scenMask));
                
                %Now preallocate a matrix in each active scenario using the
                %scenmask
                if this.calcDoseDirect
                    dij.(names{n})(this.multScen.scenMask) = {zeros(dij.doseGrid.numOfVoxels,this.numOfColumnsDij)};
                else
                    %We preallocate a sparse matrix with sparsity of
                    %1e-3 to make the filling slightly faster
                    %TODO: the preallocation could probably
                    %have more accurate estimates
                    dij.(names{n})(this.multScen.scenMask) = {spalloc(dij.doseGrid.numOfVoxels,this.numOfColumnsDij,round(prod(dij.doseGrid.numOfVoxels,this.numOfColumnsDij)*1e-3))};
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
            this.sourceModel            = this.AvailableSourceModels{1};
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
            this.dijFormatVersion       = this.defaultDijFormatVersion;

        end

        function writeHlut(this,hLutFile)

            matRad_cfg = MatRad_Config.instance();
            fileName = fullfile(this.regionsFolder, 'hLut.inp');

            mainFolder        = fullfile(matRad_cfg.matRadSrcRoot,'hluts');
            userDefinedFolder = fullfile(matRad_cfg.primaryUserFolder, 'hluts');
            fredDefinedFolder = fullfile(matRad_cfg.matRadSrcRoot, 'doseCalc', 'FRED', 'hluts');

            % Collect all the subfolders

            searchPath = [strsplit(genpath(mainFolder), pathsep)';...
                          strsplit(genpath(userDefinedFolder), pathsep)';...
                          strsplit(genpath(fredDefinedFolder),pathsep)'];
 
            searchPath(cellfun(@isempty, searchPath)) = [];

            % Check for existence of folder paths
            searchPath = searchPath(cellfun(@isfolder, searchPath));
            
            availableHLUTs = cellfun(@(x) dir([x, filesep, '*.txt']), searchPath, 'UniformOutput',false);
            availableHLUTs = cell2mat(availableHLUTs);

            hLUTindex = find(strcmp([hLutFile,'.txt'], {availableHLUTs.name}));

            if ~isempty(hLUTindex)
                selectedHlutfile = fullfile(availableHLUTs(hLUTindex).folder, availableHLUTs(hLUTindex).name);
                
                template = fileread(selectedHlutfile);
        
                newLut = fopen(fileName, 'w');
                fprintf(newLut, template);
                fclose(newLut);
            else
 
                errString = sprintf('Cannot open hLut: %s. Available hLut files are: ',hLutFile);
                errString = [errString, sprintf('\ninternal')];
                for hLUTindex=1:numel(availableHLUTs)
                    errString = [errString, sprintf('\n%s', strrep(fullfile(availableHLUTs(hLUTindex).folder, availableHLUTs(hLUTindex).name), '\', '\\'))];
                end
                matRad_cfg.dispError(errString);
            end

        end

    end

    methods (Static)
        function cmdString = cmdCall(newCmdString)
            persistent fredCmdCall;
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
                nLidx = regexp(cmdOut, '\n')+6; %6 because of tab
                nVersions = numel(nLidx)-1;

                for versIdx=1:nVersions
                    availableVersions = [availableVersions,{cmdOut(nLidx(versIdx):nLidx(versIdx)+5)}];
                end

            else
                matRad_cfg.dispError('Something wrong occured in checking FRED available version. Please check correct FRED installation');
            end
        end

        function [available,msg] = isAvailable(pln,machine)   
            % see superclass for information
            
            msg = [];
            available = false;

            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            %checkBasic
            try
                checkBasic = isfield(machine,'meta') && isfield(machine,'data');

                %check modality
                checkModality = any(strcmp(DoseEngines.matRad_ParticleFREDEngine.possibleRadiationModes, machine.meta.radiationMode));
                
                preCheck = checkBasic && checkModality;

                if ~preCheck
                    return;
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return;
            end
            available = preCheck;
        end

        function execCheck = checkExec()

            matRad_cfg = MatRad_Config.instance();
            
            %Check if I can obtain FRED version
            try
                ver = DoseEngines.matRad_ParticleFREDEngine.getVersion();
                if ~isempty(ver)
                    execCheck = true;
                else
                    execCheck = false;
                    msg = sprintf('Couldn''t call FRED executable. Please set the correct system call with DoseEngines.matRad_ParticleFREDEngine.cmdCall(''path/to/executable''). Current value is ''%s''',DoseEngines.matRad_ParticleFREDEngine.cmdCall);
                    matRad_cfg.dispError(msg);
                end
            catch
                execCheck = false;
                msg = sprintf('Couldn''t call FRED executable. Please set the correct system call with DoseEngines.matRad_ParticleFREDEngine.cmdCall(''path/to/executable''). Current value is ''%s''',DoseEngines.matRad_ParticleFREDEngine.cmdCall);
                matRad_cfg.dispError(msg);
            end
        end
        
        function version = getVersion()
            % Function to get current default FRED version
            matRad_cfg = MatRad_Config.instance();

            try
                [status, cmdOut] = system([DoseEngines.matRad_ParticleFREDEngine.cmdCall,' -vn']);

                if status == 0
                    version = cmdOut(1:end-1);
                else
                    matRad_cfg.dispError('Something wrong occured in checking FRED installation. Please check correct FRED installation');
                end
            catch
                matRad_cfg.dispWarning('Something wrong occured in checking FRED installation. Please check correct FRED installation');
                version = [];
            end
        end

        % end

        function dijMatrix = readSparseDijBin(fName)
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
            
            f = fopen(fName,'r','l');

            try
                %Header
                fileFormatVerison = fread(f,1,"int32");
                dims = fread(f,3,"int32");
                res = fread(f,3,"float32");
                offset = fread(f,3,"float32");
                nComponents = fread(f,1,"int32");
                numberOfBixels = fread(f,1,"int32");
            
                values = [];
                valuesDen = [];
                voxelIndices = [];
                colIndices = [];
                valuesNom = [];
                
                matRad_cfg.dispInfo("Reading %d number of beamlets in %d voxels (%dx%dx%d)\n",numberOfBixels,prod(dims),dims(1),dims(2),dims(3));
            
                bixelCounter = 0;
                for i = 1:numberOfBixels
                    %Read Beamlet
                    bixNum = fread(f,1,"int32");
                    numVox  = fread(f,1,"int32");
                    
                    bixelCounter = bixelCounter +1;
   
                    colIndices(end+1:end+numVox) = bixelCounter;
                    currVoxelIndices = fread(f,numVox,"uint32") + 1;
                    tmpValues = fread(f,numVox*nComponents,"float32");
                    valuesNom = tmpValues(1:nComponents:end);
                   
                    if nComponents == 2
                        valuesDen = tmpValues(nComponents:nComponents:end);
                        values(end+1:end+numVox) = valuesNom./valuesDen;
                    else
                        values(end+1:end+numVox) = valuesNom;
                    end
    
                    % x and y components have been permuted in CT
                    [indY, indX, indZ] = ind2sub(dims, currVoxelIndices);
    
                    voxelIndices(end+1:end+numVox) = sub2ind(dims([2,1,3]), indX, indY, indZ);
                    matRad_cfg.dispInfo("\tRead beamlet %d, %d voxels...\n",bixNum,numVox);
                end
                dijMatrix = sparse(voxelIndices,colIndices,values,prod(dims),numberOfBixels);
                        
                fclose(f);
            catch
                fclose(f);
                matRad_cfg.dispError('unable to load file: %s',fName);
            end
        end

        function dijMatrix = readSparseDijBin_v21(fName)
            
            matRad_cfg = MatRad_Config.instance();
            
            f = fopen(fName,'r','l');

            try
                %Header
                fileFormatVerison = fread(f,1,"int32");
                dims = fread(f,3,"int32");
                res = fread(f,3,"float32");
                offset = fread(f,3,"float32");
                orientation = fread(f,9,"float32");
                nComponents = fread(f,1,"int32");
                numberOfBixels = fread(f,1,"int32");
            
                values = [];
                valuesDen = [];
                voxelIndices = [];
                colIndices = [];
                valuesNom = [];
                
                matRad_cfg.dispInfo("Reading %d number of beamlets in %d voxels (%dx%dx%d)\n",numberOfBixels,prod(dims),dims(1),dims(2),dims(3));
            
                bixelCounter = 0;
                for i = 1:numberOfBixels
                    %Read Beamlet
                    bixNum = fread(f,1,"uint32");
                    numVox  = fread(f,1,"int32");
                    
                    bixelCounter = bixelCounter +1;
    
                    colIndices(end+1:end+numVox) = bixelCounter;
                    currVoxelIndices = fread(f,numVox,"uint32") + 1;
                    tmpValues = fread(f,numVox*nComponents,"float32");
                    valuesNom = tmpValues(1:nComponents:end);
    
                    if nComponents == 2
                        valuesDen = tmpValues(nComponents:nComponents:end);
                        values(end+1:end+numVox) = valuesNom./valuesDen;
                    else
                        values(end+1:end+numVox) = valuesNom;
                    end
    
                    % x and y components have been permuted in CT
                    [indY, indX, indZ] = ind2sub(dims, currVoxelIndices);
    
                    voxelIndices(end+1:end+numVox) = sub2ind(dims([2,1,3]), indX, indY, indZ);
                    matRad_cfg.dispInfo("\tRead beamlet %d, %d voxels...\n",bixNum,numVox);
                end
                dijMatrix = sparse(voxelIndices,colIndices,values,prod(dims),numberOfBixels);
                        
                fclose(f);
            catch
                fclose(f);
                matRad_cfg.dispError('unable to load file: %s',fName);
            end
        end

        function dijMatrix = readSparseDijBin_v31(fName)

            matRad_cfg = MatRad_Config.instance();

            f = fopen(fName,'r','l');

            try
                fileFormatVerison = fread(f,1,"int32");
                dims = fread(f,3,"int32");
                res = fread(f,3,"float32");
                offset = fread(f,3,"float32");
                orientation = fread(f,9,"float32");
                nComponents = fread(f,1,"uint32");
                numberOfBixels = fread(f,1,"uint32");
            
                values = [];
                valuesDen = [];
                voxelIndices = [];
                colIndices = [];
                valuesNom = [];
                
                matRad_cfg.dispInfo("Reading %d number of beamlets in %d voxels (%dx%dx%d)\n",numberOfBixels,prod(dims),dims(1),dims(2),dims(3));

                allBixelMeta = fread(f, 3*numberOfBixels, "uint32");
                allBixelMeta = reshape(allBixelMeta, 3,numberOfBixels)'; % (#bixels, PBidx, FID, PBID)
                
                % if nComponents>1
                %     matRad_cfg.dispWarning('!! Only last component will be read!!')
                % end
                % Components header
                for i = 1:nComponents
                    componentDataSize(i) = fread(f,1,"uint32");
                end

               % Data
               for compIdx=1:nComponents
                   % PBidx and voxelIndices should be always the same for
                   % each component?
                   PBidxs           = fread(f, componentDataSize(compIdx), "uint32")+1;
                   voxelIndices     = fread(f, componentDataSize(compIdx), "uint32")+1;
                   tmpValues{compIdx}  = fread(f, componentDataSize(compIdx), "float32");

               end
               
               % For now we only have Dose and LET scorers, if nComponents
               % == 2, then it's LET

               if nComponents>1
                   values = tmpValues{1}./tmpValues{2};
               else
                   values = tmpValues{1};
               end

               % x and y components have been permuted in CT
               [indY, indX, indZ] = ind2sub(dims, voxelIndices);

               voxelIndices = sub2ind(dims([2,1,3]), indX, indY, indZ); % + (nComponents-1)*prod(dims);
               % This will probably not work for multiple components
               dijMatrix = sparse(voxelIndices,PBidxs,values,prod(dims), numberOfBixels);

               fclose(f);
            
            catch
                fclose(f);
                matRad_cfg.dispError('unable to load file: %s',fName);
            
            end
        end

        [doseCubeV, letdCubeV, fileName] = readSimulationOutput(runFolder,calcDoseDirect, varargin);

    end


     methods (Access = private)

        function updatePaths(obj, rootFolder)

            if ~strcmp(rootFolder, obj.FREDrootFolder)
                obj.FREDrootFolder  = rootFolder;
            end
            
            obj.MCrunFolder     = fullfile(obj.FREDrootFolder, 'MCrun');
            obj.inputFolder     = fullfile(obj.MCrunFolder, 'inp');
            obj.regionsFolder   = fullfile(obj.inputFolder, 'regions');
            obj.planFolder      = fullfile(obj.inputFolder, 'plan');

        end

        function [radiationMode] = updateRadiationMode(this,value)
            % This function also resets the values for primary mass and numebr
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

        function isHigher = isVersionHigher(this,version)
            isHigher = false;
            
            % This function directly looks at FRED installation, not at
            % the current FRED version stored in the class property.
            fredVersion = this.getVersion();

            if ~isempty(fredVersion)
                % Decompose the current version for comparison
                v1 = sscanf(fredVersion, '%d.%d.%d')';
                v2 = sscanf(version, '%d.%d.%d')';
            
                if (v1(1) >= v2(1)) && (v1(2) >= v2(2)) && (v1(3) > v2(3))
                    isHigher = true;
                end
            end            

        end
     
     end

     
    methods

        function set.sourceModel(this, value)
            matRad_cfg = MatRad_Config.instance();

            valid = ischar(value) && any(strcmp(value, this.AvailableSourceModels));

            if valid
                this.sourceModel = value;
            else
                matRad_cfg.dispWarning('Unable to set source model:%s, setting default:%s', value, this.AvailableSourceModels{1})
                this.sourceModel = this.AvailableSourceModels{1};
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

         function set.dijFormatVersion(this,value)

             matRad_cfg = MatRad_Config.instance();

             if ~this.isVersionHigher('3.70.0')
                % FRED version < 3.70.0 does not allow dij version
                % selection and only works with ifFormatVersion < 21
   
                this.dijFormatVersion = this.defaultDijFormatVersion;
             
                if ~strcmp(value, this.defaultDijFormatVersion)
                    matRad_cfg.dispWarning(sprintf('ijFormat: %s not available for FRED version < 3.70.0', value));
                end
             
             else
                 %if this.useGPU
                     if strcmp(value, '20')
                         this.dijFormatVersion = '21';
                     else
                         this.dijFormatVersion = value;
                     end
                 %else
                 %    this.dijFormatVersion = '20';
                 %end
             end
         end

         function readerHandle = get.dijReaderHandle(this)

             matRad_cfg = MatRad_Config.instance();

             switch this.dijFormatVersion

                 case '20'
                    readerHandle = @(lFile) DoseEngines.matRad_ParticleFREDEngine.readSparseDijBin(lFile);                     
                 case '21'
                    readerHandle = @(lFile) DoseEngines.matRad_ParticleFREDEngine.readSparseDijBin_v21(lFile);
                 case '31'
                    readerHandle = @(lFile) DoseEngines.matRad_ParticleFREDEngine.readSparseDijBin_v31(lFile);
                 otherwise
                     matRad_cfg.dispWarning(sprintf('Unable to read dij format version: %s, using default: %s', this.dijFormatVersion, this.defaultDijFormatVersion));
                     readerHandle = @(lFile) DoseEngines.matRad_ParticleFREDEngine.readSparseDijBin(lFile);
             end

                
         end

         function set.radiationMode(this,value)
             if ischar(value)
                 if ~isempty(this.radiationMode) && ~strcmp(this.radiationMode, value)
                    this.radiationMode = value;
                    this.updateRadiationMode(this.radiationMode);
                 elseif isempty(this.radiationMode)
                    this.radiationMode = value;
                 end
             end

         end

        function set.FREDrootFolder(obj, pathValue)
            obj.FREDrootFolder = pathValue;
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

