classdef matRad_ParticleMCsquareEngine < DoseEngines.matRad_MonteCarloEngineAbstract
% Engine for particle dose calculation using monte carlo calculation
% specificly the mc square method
% for more informations see superclass
% DoseEngines.matRad_MonteCarloEngineAbstract
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    properties (Constant)  
        possibleRadiationModes = {'protons'};
        name = 'MCsquare';
        shortName = 'MCsquare';
    end

    properties
        config;             %Holds an instance of all configurable parameters (matRad_MCsquareConfig)
        MCsquareFolder;     %Folder to the MCsquare installation
        forceBDL = [];      %Specify an existing BDL file to load

        %Other Dose Calculation Properties
        calcLET = false;        
    end
    
    properties (SetAccess = protected, GetAccess = public)
        
        currFolder = pwd; %folder path when set
        
        mcSquareBinary; %Executable for mcSquare simulation
        nbThreads; %number of threads for MCsquare, 0 is all available

        constantRBE = NaN;              % constant RBE value 
    end      
    
    methods
        
        function this = matRad_ParticleMCsquareEngine(pln)
            % Constructor
            %
            % call
            %   engine = DoseEngines.matRad_DoseEngineMCsquare(ct,stf,pln,cst)
            %
            % input
            %   ct:                         matRad ct struct
            %   stf:                        matRad steering information struct
            %   pln:                        matRad plan meta information struct
            %   cst:                        matRad cst struct
            
            if nargin < 1
                pln = [];
            end

            % call superclass constructor
            this = this@DoseEngines.matRad_MonteCarloEngineAbstract(pln);

            this.config = matRad_MCsquareConfig();
            
            % check if bio optimization is needed and set the
            % coresponding boolean accordingly
            % TODO:
            % This should not be handled here as an optimization property
            % We should rather make optimization dependent on what we have
            % decided to calculate here.
            if nargin > 0 
                if (isfield(pln,'propOpt')&& isfield(pln.propOpt,'bioOptimization')&& ...
                    (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') ||...
                    isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) && ...
                    strcmp(pln.radiationMode,'carbon'))
                this.calcBioDose = true;
                elseif strcmp(pln.radiationMode,'protons') && isfield(pln,'propOpt') && isfield(pln.propOpt,'bioOptimization') && isequal(pln.propOpt.bioOptimization,'const_RBExD')
                    this.constantRBE = 1.1;                    
                end
            end
        end

        function setDefaults(this)
            this.setDefaults@DoseEngines.matRad_MonteCarloEngineAbstract();

            % future code for property validation on creation here
            matRad_cfg = MatRad_Config.instance();
            
            %Assign default parameters from MatRad_Config
            this.doseGrid.resolution    = matRad_cfg.propDoseCalc.defaultResolution;
            this.multScen = 'nomScen';
            this.selectVoxelsInScenarios = matRad_cfg.propDoseCalc.defaultSelectVoxelsInScenarios;

            %Set Default MCsquare path
            %Set folder
            this.MCsquareFolder = [matRad_cfg.matRadRoot filesep 'MCsquare' filesep 'bin'];
        end
    end
    
    methods(Access = protected)
        
        function dij = calcDose(this,ct,cst,stf)
            % matRad MCsqaure monte carlo photon dose calculation wrapper
            % can be automaticly called through matRad_calcDose or
            % matRad_calcParticleDoseMC
            %
            % nCase per Bixel and be either set by hand after creating the
            % engine or over the matRad_calcPhotonDoseMC function while
            % calling the calculation
            %
            % call
            %   dij = this.calcDose(ct,stf,pln,cst)
            %
            % input
            %   ct:          	matRad ct struct
            %   cst:            matRad cst struct
            %   stf:         	atRad steering information struct
            %   
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

            %Now we can run initDoseCalc as usual
            dij = this.initDoseCalc(ct,cst,stf);

            % switch for using existing BDL file (e.g. to fit matRad basedata),
            % or generate BDL file from matRad base data using MCsquareBDL
            if ~isempty(this.forceBDL)
                % use existing BDL file
                bdFile = this.forceBDL;

            else
                % fit and create BDL file using selected machine file
                bdFile = [this.machine.meta.machine '.txt'];

                % Calculate MCsquare base data
                % Argument stf is optional, if given, calculation only for energies given in stf
                MCsquareBDL = matRad_MCsquareBaseData(this.machine);

                %matRad_createMCsquareBaseDataFile(bdFile,machine,1);
                MCsquareBDL = MCsquareBDL.writeMCsquareData([this.MCsquareFolder filesep 'BDL' filesep bdFile]);
                MCsquareBDL = MCsquareBDL.saveMatradMachine('savedMatRadMachine');

            end

            
            % The offset of the dose grid of MCsquare and matRad are
            % different, so we can not use the one computed in
            % dij.doseGrid.isoCenterOffset
            mcSquareAddIsoCenterOffset = [dij.doseGrid.resolution.x/2 dij.doseGrid.resolution.y/2 dij.doseGrid.resolution.z/2] ...
                            - [dij.ctGrid.resolution.x   dij.ctGrid.resolution.y   dij.ctGrid.resolution.z];

            % for MCsquare we explicitly downsample the ct to the dose grid (might not
            % be necessary in future MCsquare versions with separated grids)
            for s = 1:dij.numOfScenarios
                HUcube{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cubeHU{s}, ...
                                            dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
            end
           
            % set absolute calibration factor
            % convert from eV/g/primary to Gy 1e6 primaries
            absCalibrationFactorMC2 = 1.602176e-19 * 1.0e+9;


            % MCsquare settings
            MCsquareConfigFile = 'MCsquareConfig.txt';

            this.config.BDL_Plan_File = 'currBixels.txt'; 
            this.config.BDL_Machine_Parameter_File = ['BDL/' bdFile];

            this.config.CT_File       = 'MC2patientCT.mhd';
            this.config.Num_Threads   = this.nbThreads;
            this.config.RNG_Seed      = 1234;
            if this.calcDoseDirect 
                this.config.Num_Primaries = this.numHistoriesDirect;
            else
                this.config.Num_Primaries = this.numHistoriesPerBeamlet;
            end

            % turn simulation of individual beamlets
            this.config.Beamlet_Mode = ~this.calcDoseDirect;
            % turn of writing of full dose cube
            this.config.Dose_MHD_Output = this.calcDoseDirect;
            % turn on sparse output
            this.config.Dose_Sparse_Output = ~this.calcDoseDirect;
            % set threshold of sparse matrix generation
            this.config.Dose_Sparse_Threshold = 1 - this.relativeDosimetricCutOff;

            %Matrices for LET
            if this.calcLET
                this.config.LET_MHD_Output		 = this.calcDoseDirect;
                this.config.LET_Sparse_Output	 = ~this.calcDoseDirect;
            end

            for scenarioIx = 1:this.multScen.totNumScen
                %For direct dose calculation
                totalWeights = 0;

                % manipulate isocenter
                isoCenterShift = this.multScen.isoShift(scenarioIx,:) + mcSquareAddIsoCenterOffset;
                
                ctScen = this.multScen.linearMask(scenarioIx,1);
                shiftScen = this.multScen.linearMask(scenarioIx,2);
                rangeShiftScen = this.multScen.linearMask(scenarioIx,3);

                if this.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
           
 
                    counter = 0;
                    for i = 1:length(stf)
                        %Create new stf for MCsquare with energy layer ordering and
                        %shifted scenario isocenter
                        stfMCsquare(i).isoCenter   = stf(i).isoCenter + isoCenterShift;
                        stfMCsquare(i).gantryAngle = mod(180-stf(i).gantryAngle,360); %Different MCsquare geometry
                        stfMCsquare(i).couchAngle  = stf(i).couchAngle;
                        stfMCsquare(i).energies    = unique([stf(i).ray.energy]);
                        stfMCsquare(i).SAD         = stf(i).SAD;

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

                        % allocate empty target point container
                        for j = 1:numel(stfMCsquare(i).energies)
                            stfMCsquare(i).energyLayer(j).targetPoints   = [];
                            stfMCsquare(i).energyLayer(j).numOfPrimaries = [];
                            stfMCsquare(i).energyLayer(j).MU             = [];
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
                                    energyIx = find(stf(i).ray(j).energy == stfMCsquare(i).energies(k));
                                    stfMCsquare(i).energyLayer(k).rayNum   = [stfMCsquare(i).energyLayer(k).rayNum j];
                                    stfMCsquare(i).energyLayer(k).bixelNum = [stfMCsquare(i).energyLayer(k).bixelNum energyIx];
                                    stfMCsquare(i).energyLayer(k).targetPoints = [stfMCsquare(i).energyLayer(k).targetPoints; ...
                                        -stf(i).ray(j).rayPos_bev(1) stf(i).ray(j).rayPos_bev(3)];

                                    %Number of primaries depending on beamlet-wise or field-based compuation (direct dose calculation)
                                    if this.calcDoseDirect
                                        stfMCsquare(i).energyLayer(k).numOfPrimaries = [stfMCsquare(i).energyLayer(k).numOfPrimaries ...
                                            round(stf(i).ray(j).weight(stf(i).ray(j).energy == stfMCsquare(i).energies(k))*this.numHistoriesDirect)];

                                        stfMCsquare(i).energyLayer(k).MU = [stfMCsquare(i).energyLayer(k).MU ...
                                            round(stf(i).ray(j).weight(stf(i).ray(j).energy == stfMCsquare(i).energies(k))*this.numHistoriesDirect)];

                                        totalWeights = totalWeights + stf(i).ray(j).weight(stf(i).ray(j).energy == stfMCsquare(i).energies(k));
                                    else
                                        stfMCsquare(i).energyLayer(k).numOfPrimaries = [stfMCsquare(i).energyLayer(k).numOfPrimaries ...
                                            this.numHistoriesPerBeamlet];

                                        stfMCsquare(i).energyLayer(k).MU = [stfMCsquare(i).energyLayer(k).MU ...
                                            this.numHistoriesPerBeamlet];
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
                    % write patient data
                    MCsquareBinCubeResolution = [dij.doseGrid.resolution.x ...
                        dij.doseGrid.resolution.y ...
                        dij.doseGrid.resolution.z];

                    this.writeMhd(HUcube{ctScen},MCsquareBinCubeResolution);

                    % write config file
                    this.writeInputFiles(MCsquareConfigFile,stfMCsquare);

                    %% MC computation and dij filling
                    

                    % run MCsquare
                    mcSquareCall = [this.mcSquareBinary ' ' MCsquareConfigFile];
                    matRad_cfg.dispInfo(['Calling Monte Carlo Engine: ' mcSquareCall]);
                    [status,cmdout] = system(mcSquareCall,'-echo');

                    mask = false(dij.doseGrid.numOfVoxels,1);
                    mask(this.VdoseGrid) = true;

                    % read sparse matrix
                    if ~this.calcDoseDirect
                        dij.physicalDose{ctScen,shiftScen,rangeShiftScen} = absCalibrationFactorMC2 * matRad_sparseBeamletsReaderMCsquare ( ...
                            [this.config.Output_Directory filesep 'Sparse_Dose.bin'], ...
                            dij.doseGrid.dimensions, ...
                            dij.totalNumOfBixels, ...
                            mask);
                        
                        %Read sparse LET
                        if this.calcLET
                            dij.mLETDose{ctScen,shiftScen,rangeShiftScen} = dij.physicalDose{ctScen,shiftScen,rangeShiftScen} .* matRad_sparseBeamletsReaderMCsquare ( ...
                                [this.config.Output_Directory filesep 'Sparse_LET.bin'], ...
                                dij.doseGrid.dimensions, ...
                                dij.totalNumOfBixels, ...
                                mask);
                        end
                        
                        % reorder influence matrix to comply with matRad default ordering
                        dij.physicalDose = cellfun(@(mx) mx(:,MCsquareOrder),dij.physicalDose,'UniformOutput',false);
                        if this.calcLET
                            dij.mLETDose = cellfun(@(mx) mx(:,MCsquareOrder),dij.mLETDose,'UniformOutput',false);
                        end
                    else
                        cube = this.readMhd('Dose.mhd');
                        dij.physicalDose{ctScen,shiftScen,rangeShiftScen} = sparse(this.VdoseGrid,ones(numel(this.VdoseGrid),1), ...
                            absCalibrationFactorMC2 * cube(this.VdoseGrid), ...
                            dij.doseGrid.numOfVoxels,1);

                        %Read LET cube
                        if this.calcLET
                            cube = this.readMhd('LET.mhd');
                            dij.mLETDose{ctScen,shiftScen,rangeShiftScen} = dij.physicalDose{ctScen,shiftScen,rangeShiftScen} .* sparse(this.VdoseGrid,ones(numel(this.VdoseGrid),1), ...
                                cube(this.VdoseGrid), ...
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

                    
                    if this.config.Beamlet_Mode
                        
                    end

                    matRad_cfg.dispInfo('Scenario %d of %d finished!\n',scenarioIx,this.multScen.totNumScen);
                    
                    %% clear all data
                    %could also be moved to the "finalize" function
                    delete([this.config.CT_File(1:end-4) '.*']);
                    delete('currBixels.txt');
                    delete('MCsquareConfig.txt');

                    %For Octave temporarily disable confirmation for recursive rmdir
                    if strcmp(matRad_cfg.env,'OCTAVE')
                        rmdirConfirmState = confirm_recursive_rmdir(0);
                    end
                    rmdir(this.config.Output_Directory,'s');

                    %Reset to old confirmatoin state
                    if strcmp(matRad_cfg.env,'OCTAVE')
                        confirm_recursive_rmdir(rmdirConfirmState);
                    end

                    % cd back
                    cd(this.currFolder);
                end
            end
        
            matRad_cfg.dispInfo('matRad: Simulation finished!\n');
            %Finalize dose calculation
            dij = this.finalizeDose(dij);

        end

        function setBinaries(this)
            % setBinaries check if the binaries are available on the current
            % machine and sets to the mcsquarebinary object property
            %

            [~,binaryFile] = this.checkBinaries();
            this.mcSquareBinary = binaryFile;
        end
        
        function dij = initDoseCalc(this,ct,cst,stf)
            %% Assingn and check parameters
            matRad_cfg = MatRad_Config.instance();            

            % check if binaries are available
            % Executables for simulation
            this.setBinaries();

            %Mex interface for import of sparse matrix
            if ~this.calcDoseDirect && ~matRad_checkMexFileExists('matRad_sparseBeamletsReaderMCsquare')
                matRad_cfg.dispWarning('Compiled sparse reader interface not found. Trying to compile it on the fly!');      
                try
                    matRad_compileMCsquareSparseReader();        
                catch MException
                    matRad_cfg.dispError('Could not find/generate mex interface for reading the sparse matrix. \nCause of error:\n%s\n Please compile it yourself.',MException.message);
                end
            end

            % set and change to MCsquare binary folder
            this.currFolder = pwd;
            %fullfilename = mfilename('fullpath');

            % cd to MCsquare folder (necessary for binary)
            % TODO: Could be checked in a property setter function
            if ~exist(this.MCsquareFolder,'dir')
                matRad_cfg.dispError('MCsquare Folder does not exist!'); 
            end
            cd(this.MCsquareFolder);

            %Check Materials
            if ~exist([this.MCsquareFolder filesep 'Materials'],'dir') || ~exist(fullfile(this.MCsquareFolder,'Materials','list.dat'),'file')
                matRad_cfg.dispInfo('First call of MCsquare: unzipping Materials...');    
                unzip('Materials.zip');
                matRad_cfg.dispInfo('Done');
            end

            %% Call Superclass init function
            dij = initDoseCalc@DoseEngines.matRad_MonteCarloEngineAbstract(this,ct,cst,stf); 

            % Since MCsquare 1.1 only allows similar resolution in x&y, we do some
            % extra checks on that before calling the normal initDoseCalc. First, we make sure a
            % dose grid resolution is set in the pln struct

            % Now we check for different x/y
            if dij.doseGrid.resolution.x ~= dij.doseGrid.resolution.y
                dij.doseGrid.resolution.x = mean([dij.doseGrid.resolution.x dij.doseGrid.resolution.y]);
                dij.doseGrid.resolution.y = dij.doseGrid.resolution.x;
                matRad_cfg.dispWarning('Anisotropic resolution in axial plane for dose calculation with MCsquare not possible\nUsing average x = y = %g mm\n',dij.doseGrid.resolution.x);
            end
            
            %% Validate and preset some additional dij variables
            
            % Explicitly setting the number of threads for MCsquare, 0 is all available
            this.nbThreads = 0;
            
            %Issue a warning when we have more than 1 scenario
            if dij.numOfScenarios ~= 1
                matRad_cfg.dispWarning('MCsquare is only implemented for single scenario use at the moment. Will only use the first Scenario for Monte Carlo calculation!');
            end
            
            % prefill ordering of MCsquare bixels
            dij.MCsquareCalcOrder = NaN*ones(dij.totalNumOfBixels,1);  

            if ~isnan(this.constantRBE) 
                dij.RBE = this.constantRBE;
            end
        end
        
        function writeInputFiles(obj,filename,stf)           
            % generate input files for MCsquare dose calcualtion from matRad
            %
            % call
            %   obj.writeInputFiles(filename,filename,stf)
            %
            % input
            %   filename:       filename of the Configuration file
            %   stf:            matRad steering information struct
            %
            % output
            %   -
            %
            % References
            %   [1] https://openreggui.org/git/open/REGGUI/blob/master/functions/io/convert_Plan_PBS.m
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
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %% write overall configuration file
            fileHandle = fopen(filename,'w');
            obj.config.write(fileHandle);
            fclose(fileHandle);

            %% prepare steering file writing
            numOfFields = length(stf);
            if obj.config.Beamlet_Mode
                totalMetersetWeightOfAllFields = 1;
            else
                totalMetersetWeightOfFields = NaN*ones(numOfFields,1);
                for i = 1:numOfFields
                    totalMetersetWeightOfFields(i) = sum([stf(i).energyLayer.numOfPrimaries]);
                end
                totalMetersetWeightOfAllFields = sum(totalMetersetWeightOfFields);
            end

            %% write steering file

            fileHandle = fopen(obj.config.BDL_Plan_File,'w');

            fprintf(fileHandle,'#TREATMENT-PLAN-DESCRIPTION\n');
            fprintf(fileHandle,'#PlanName\n');
            fprintf(fileHandle,'matRad_bixel\n');
            fprintf(fileHandle,'#NumberOfFractions\n');
            fprintf(fileHandle,'1\n');
            fprintf(fileHandle,'##FractionID\n');
            fprintf(fileHandle,'1\n');
            fprintf(fileHandle,'##NumberOfFields\n');
            fprintf(fileHandle,[num2str(numOfFields) '\n']);
            for i = 1:numOfFields
                fprintf(fileHandle,'###FieldsID\n');
                fprintf(fileHandle,[num2str(i) '\n']);
            end
            fprintf(fileHandle,'\n#TotalMetersetWeightOfAllFields\n');
            fprintf(fileHandle,[num2str(totalMetersetWeightOfAllFields) '\n']);

            for i = 1:numOfFields
                fprintf(fileHandle,'\n#FIELD-DESCRIPTION\n');
                fprintf(fileHandle,'###FieldID\n');
                fprintf(fileHandle,[num2str(i) '\n']);
                fprintf(fileHandle,'###FinalCumulativeMeterSetWeight\n');
                if obj.config.Beamlet_Mode
                    finalCumulativeMeterSetWeight = 1/numOfFields;
                else
                    finalCumulativeMeterSetWeight = totalMetersetWeightOfFields(i);
                end
                fprintf(fileHandle,[num2str(finalCumulativeMeterSetWeight) '\n']);
                fprintf(fileHandle,'###GantryAngle\n');
                fprintf(fileHandle,[num2str(stf(i).gantryAngle) '\n']);
                fprintf(fileHandle,'###PatientSupportAngle\n');
                fprintf(fileHandle,[num2str(stf(i).couchAngle) '\n']);
                fprintf(fileHandle,'###IsocenterPosition\n');
                fprintf(fileHandle,[num2str(stf(i).isoCenter) '\n']);
                fprintf(fileHandle,'###NumberOfControlPoints\n');
                numOfEnergies = numel(stf(i).energies);
                fprintf(fileHandle,[num2str(numOfEnergies) '\n']);

                %Range shfiter
                if stf(i).rangeShifterID ~= 0
                    fprintf(fileHandle,'###RangeShifterID\n%d\n',stf(i).rangeShifterID);
                    fprintf(fileHandle,'###RangeShifterType\n%s\n',stf(i).rangeShifterType);
                end

                metersetOffset = 0;
                fprintf(fileHandle,'\n#SPOTS-DESCRIPTION\n');
                for j = 1:numOfEnergies
                    fprintf(fileHandle,'####ControlPointIndex\n');
                    fprintf(fileHandle,[num2str(j) '\n']);
                    fprintf(fileHandle,'####SpotTunnedID\n');
                    fprintf(fileHandle,['1\n']);
                    fprintf(fileHandle,'####CumulativeMetersetWeight\n');
                    if obj.config.Beamlet_Mode
                        cumulativeMetersetWeight = j/numOfEnergies * 1/numOfFields;
                    else
                        cumulativeMetersetWeight = metersetOffset + sum([stf(i).energyLayer(j).numOfPrimaries]);
                        metersetOffset = cumulativeMetersetWeight;
                    end
                    fprintf(fileHandle,[num2str(cumulativeMetersetWeight) '\n']);
                    fprintf(fileHandle,'####Energy (MeV)\n');
                    fprintf(fileHandle,[num2str(stf(i).energies(j)) '\n']);

                    %Range shfiter
                    if stf(i).rangeShifterID ~= 0
                        rangeShifter = stf(i).energyLayer(j).rangeShifter;
                        if rangeShifter.ID ~= 0
                            fprintf(fileHandle,'####RangeShifterSetting\n%s\n','IN');
                            pmma_rsp = 1.165; %TODO: hardcoded for now
                            rsWidth = rangeShifter.eqThickness / pmma_rsp;
                            isoToRaShi = stf(i).SAD - rangeShifter.sourceRashiDistance + rsWidth;
                            fprintf(fileHandle,'####IsocenterToRangeShifterDistance\n%f\n',-isoToRaShi/10); %in cm
                            fprintf(fileHandle,'####RangeShifterWaterEquivalentThickness\n%f\n',rangeShifter.eqThickness);
                        else
                            fprintf(fileHandle,'####RangeShifterSetting\n%s\n','OUT');
                        end
                    end

                    fprintf(fileHandle,'####NbOfScannedSpots\n');
                    numOfSpots = size(stf(i).energyLayer(j).targetPoints,1);
                    fprintf(fileHandle,[num2str(numOfSpots) '\n']);
                    fprintf(fileHandle,'####X Y Weight\n');
                    for k = 1:numOfSpots
                        %{
                        if obj.config.Beamlet_Mode
                            n = stf(i).energyLayer(j).numOfPrimaries(k);
                        else
                            n = stf(i).energyLayer(j).numOfPrimaries(k) /  obj.mcSquare_magicFudge(stf(i).energies(j));
                        end
                        %}
                        n = stf(i).energyLayer(j).numOfPrimaries(k);
                        fprintf(fileHandle,[num2str(stf(i).energyLayer(j).targetPoints(k,:)) ' ' num2str(n) '\n']);
                    end
                end
            end

            fclose(fileHandle);

        end

        function cube = readMhd(obj,filename)
            % TODO: This should become a binary export function in matRads
            % IO folde
            % matRad mhd file reader
            %
            % call
            %   cube = matRad_readMhd(folder,filename)
            %
            % input
            %   folder:   folder where the *raw and *mhd file are located
            %   filename: filename
            %
            % output
            %   cube:     3D array
            %
            % References
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
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %% read header
            headerFileHandle = fopen([obj.config.Output_Directory, filesep filename],'r');

            s = textscan(headerFileHandle, '%s', 'delimiter', '\n');

            % read dimensions
            idx = find(~cellfun(@isempty,strfind(s{1}, 'DimSize')),1,'first');
            dimensions = cell2mat(textscan(s{1}{idx},'DimSize = %f %f %f'));

            % read filename of data
            idx = find(~cellfun(@isempty,strfind(s{1}, 'ElementDataFile')),1,'first');
            tmp = textscan(s{1}{idx},'ElementDataFile = %s');
            dataFilename = cell2mat(tmp{1});

            % get data type
            idx = find(~cellfun(@isempty,strfind(s{1}, 'ElementType')),1,'first');
            tmp = textscan(s{1}{idx},'ElementType = MET_%s');
            type = lower(cell2mat(tmp{1}));

            fclose(headerFileHandle);

            %% read data
            dataFileHandle = fopen([obj.config.Output_Directory filesep dataFilename],'r');
            cube = reshape(fread(dataFileHandle,inf,type),dimensions);
            cube = permute(cube,[2 1 3]);
            cube = flip(cube,2);
            fclose(dataFileHandle);
        end
    
        function writeMhd(obj,cube,resolution)
            % TODO: This should become a binary export function in matRads
            % IO folder
            % References
            %   -
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Copyright 2020 the matRad development team.
            %
            % This file is part of the matRad project. It is subject to the license
            % terms in the LICENSE file found in the top-level directory of this
            % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
            % of the matRad project, including this file, may be copied, modified,
            % propagated, or distributed except according to the terms contained in the
            % LICENSE file.
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% write header file
            fileHandle = fopen(obj.config.CT_File,'w');

            fprintf(fileHandle,'ObjectType = Image\n');
            fprintf(fileHandle,'NDims = 3\n');
            fprintf(fileHandle,'BinaryData = True\n');
            fprintf(fileHandle,'BinaryDataByteOrderMSB = False\n');
            fprintf(fileHandle,'CompressedData = False\n');
            fprintf(fileHandle,'TransformMatrix = 1 0 0 0 1 0 0 0 1\n');
            fprintf(fileHandle,'Offset = 0 0 0\n');
            fprintf(fileHandle,'CenterOfRotation = 0 0 0\n');
            fprintf(fileHandle,'AnatomicalOrientation = RAI\n');
            fprintf(fileHandle,'ElementSpacing = %f %f %f\n',resolution);
            fprintf(fileHandle,'DimSize = %d %d %d\n',size(cube,2),size(cube,1),size(cube,3));
            fprintf(fileHandle,'ElementType = MET_DOUBLE\n');
            filenameRaw = [obj.config.CT_File(1:end-4) '.raw'];
            fprintf(fileHandle,'ElementDataFile = %s\n',filenameRaw);

            fclose(fileHandle);

            %% write data file
            dataFileHandle = fopen(filenameRaw,'w');

            cube = flip(cube,2);
            cube = permute(cube,[2 1 3]);

            fwrite(dataFileHandle,cube(:),'double');
            fclose(dataFileHandle);
        end

    end

    methods (Access = private)
        function gain = mcSquare_magicFudge(~,energy)
            % mcSquare will scale the spot intensities in
            % https://gitlab.com/openmcsquare/MCsquare/blob/master/src/data_beam_model.c#L906
            % by this factor so we need to divide up front to make things work. The
            % original code can be found at https://gitlab.com/openmcsquare/MCsquare/blob/master/src/compute_beam_model.c#L16

            K = 35.87; % in eV (other value 34.23 ?)

            % // Air stopping power (fit ICRU) multiplied by air density
            SP = (9.6139e-9*energy^4 - 7.0508e-6*energy^3 + 2.0028e-3*energy^2 - 2.7615e-1*energy + 2.0082e1) * 1.20479E-3 * 1E6; % // in eV / cm

            % // Temp & Pressure correction
            PTP = 1.0;

            % // MU calibration (1 MU = 3 nC/cm)
            % // 1cm de gap effectif
            C = 3.0E-9; % // in C / cm

            % // Gain: 1eV = 1.602176E-19 J
            gain = (C*K) / (SP*PTP*1.602176E-19);

            % divide by 1e7 to not get tiny numbers...
            gain = gain/1e7;

        end
    end

    
    methods (Static)

        function [binaryFound,binaryFile] = checkBinaries()
            % checkBinaries check if the binaries are available on the current
            % machine and sets to the mcsquarebinary object property
            %        
            %
            matRad_cfg = MatRad_Config.instance();
            
            binaryFile = [];
            binaryFound = false;

            if ispc
                if exist('MCSquare_windows.exe','file') ~= 2
                    matRad_cfg.dispWarning('Could not find MCsquare binary.\n');
                else
                    binaryFile = 'MCSquare_windows.exe';
                end
            elseif ismac
                if exist('MCsquare_mac','file') ~= 2
                    matRad_cfg.dispWarning('Could not find MCsquare binary.\n');
                else
                    binaryFile = './MCsquare_mac';
                end
                %error('MCsquare binaries not available for mac OS.\n');
            elseif isunix
                if exist('MCsquare_linux','file') ~= 2
                    matRad_cfg.dispWarning('Could not find MCsquare binary.\n');
                else
                    binaryFile = './MCsquare_linux';
                end
            end

            if ~isempty(binaryFile)
                binaryFound = true;
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
                checkModality = any(strcmp(DoseEngines.matRad_ParticleMCsquareEngine.possibleRadiationModes, machine.meta.radiationMode));
                
                preCheck = checkBasic && checkModality;

                if ~preCheck
                    return;
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return;
            end

            %Check the binaries
            hasBinaries = DoseEngines.matRad_ParticleMCsquareEngine.checkBinaries();

            available = preCheck & hasBinaries;
        end
              
    end
end

