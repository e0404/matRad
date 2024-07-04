classdef matRad_TopasMCEngine < DoseEngines.matRad_MonteCarloEngineAbstract
    % matRad_TopasMCEngine 
    %   Implementation of the TOPAS interface for Monte Carlo dose
    %   calculation
    %
    % References
    %
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
        possibleRadiationModes = {'photons','protons','helium','carbon'};
        name = 'TOPAS';
        shortName = 'TOPAS';

        defaultPhotonBeamProfile = 'uniform';
        defaultParticleBeamProfile = 'biGaussian';
    end

    properties
        calcLET = false;
        calcBioDose = false;
        prescribedDose = [];


        topasExecCommand; %Defaults will be set during construction according to TOPAS installation instructions and used system

        parallelRuns = false; %Starts runs in parallel
        
        externalCalculation = 'off'; %Generates folder for external TOPAS calculation (e.g. on a server)

        workingDir; %working directory for the simulation

        engine = 'TOPAS'; %parameter for continuity

        label = 'matRad_plan';

        %Simulation parameters
        numThreads = 0; %number of used threads, 0 = max number of threads (= num cores)
        numOfRuns = 1; %Default number of runs / batches
        modeHistories = 'num'; %'frac';
        fracHistories = 1e-4; %Fraction of histories to compute

        numParticlesPerHistory = 1e6;
        verbosity = struct( 'timefeatures',0,...
            'cputime',true,...
            'run',0,...
            'event',0,...
            'tracking',0,...
            'material',0,...
            'maxinterruptedhistories',1000,...
            'maxDetailedErrorReports',0);

        minRelWeight = .00001; %Threshold for discarding beamlets. 0 means all weights are being considered, can otherwise be assigned to min(w)

        useOrigBaseData = false; % base data of the original matRad plan will be used?
        beamProfile = 'biGaussian'; %'biGaussian' (emittance); 'simple'

        useEnergySpectrum = false;

        %Not yet implemented
        %beamletMode = false; %In beamlet mode simulation will be performed for a dose influence matrix (i.e., each beamlet simulates numHistories beamlets)

        pencilBeamScanning = true; %This should be always true except when using photons (enables deflection)

        %Image
        materialConverter = struct('mode','HUToWaterSchneider',...    %'RSP','HUToWaterSchneider';
            'densityCorrection','Schneider_TOPAS',... %'rspHLUT','Schneider_TOPAS','Schneider_matRad'
            'addSection','none',... %'none','lung'
            'addTitanium',false,... %'false','true' (can only be used with advanced HUsections)
            'HUSection','advanced',... %'default','advanced'
            'HUToMaterial','default',... %'default',','advanced','MCsquare'
            'loadConverterFromFile',false); % set true if you want to use your own SchneiderConverter written in "TOPAS_SchneiderConverter"

        arrayOrdering = 'F'; %'C';
        rsp_basematerial = 'Water';

        %Scoring
        scorer = struct('volume',false,...
            'doseToMedium',true,...
            'doseToWater',false,...
            'surfaceTrackCount',false,...
            'calcDij',false,...
            'RBE',false,...
            'RBE_model',{{'default'}},... % default is MCN for protons and LEM1 for ions
            'defaultModelProtons',{{'MCN'}},...
            'defaultModelCarbon',{{'LEM'}},...
            'LET',false,...
            'sharedSubscorers',true,...
            'outputType','binary',... %'csv'; 'binary';%
            ... % This variable is only used for physicalDose, since for now it adds unnecessary computation time
            'reportQuantity',{{'Sum','Standard_Deviation'}});         % 'reportQuantity',{{'Sum'}});
        scorerRBEmodelOrderForEvaluation = {'MCN','WED','LEM','libamtrack'};
        bioParameters = struct( 'PrescribedDose',2,...
            'AlphaX',0.1,...
            'BetaX',0.05,...
            'SimultaneousExposure','"True"');

        %Physics
        electronProductionCut = 0.5; %in mm
        radiationMode;
        modules_protons     = {'g4em-standard_opt4','g4h-phy_QGSP_BIC_HP','g4decay','g4h-elastic_HP','g4stopping','g4ion-QMD','g4radioactivedecay'};
        modules_GenericIon  = {'g4em-standard_opt4','g4h-phy_QGSP_BIC_HP','g4decay','g4h-elastic_HP','g4stopping','g4ion-QMD','g4radioactivedecay'};
        modules_photons       = {'g4em-standard_opt4','g4h-phy_QGSP_BIC_HP','g4decay'};

        %Geometry / World
        worldMaterial = 'G4_AIR';

        %filenames
        converterFolder = 'materialConverter';
        scorerFolder = 'scorer';
        outfilenames = struct(  'patientParam','matRad_cube.txt',...
            'patientCube','matRad_cube.dat');

        infilenames = struct(   'geometry','world/TOPAS_matRad_geometry.txt.in',...
            ... % BeamSetup files
            'beam_virtualGaussian','beamSetup/TOPAS_beamSetup_virtualGaussian.txt.in',...
            'beam_phasespace','beamSetup/TOPAS_beamSetup_phasespace.txt.in',...
            'beam_uniform','beamSetup/TOPAS_beamSetup_uniform.txt.in',...
            'beam_mlc','beamSetup/TOPAS_beamSetup_mlc.txt.in',...
            'beam_biGaussian','beamSetup/TOPAS_beamSetup_biGaussian.txt.in',...
            'beam_generic','beamSetup/TOPAS_beamSetup_generic.txt.in',...
            ... % Schneier Converter
            ... % Defined Materials
            'matConv_Schneider_definedMaterials',struct('default','definedMaterials/default.txt.in',...
            'MCsquare','definedMaterials/MCsquare.txt.in',...
            'advanced','definedMaterials/advanced.txt.in'),...
            ... % Density Correction
            'matConv_Schneider_densityCorr_Schneider_matRad','densityCorrection/Schneider_matRad.dat',...
            'matConv_Schneider_densityCorr_Schneider_TOPAS','densityCorrection/Schneider_TOPAS.dat',...
            ... % load from file
            'matConv_Schneider_loadFromFile','TOPAS_SchneiderConverter.txt.in',...
            ... % Scorer
            'Scorer_surfaceTrackCount','TOPAS_scorer_surfaceIC.txt.in',...
            'Scorer_doseToMedium','TOPAS_scorer_doseToMedium.txt.in',...
            'Scorer_LET','TOPAS_subscorer_LET.txt.in',...
            'Scorer_doseToWater','TOPAS_scorer_doseToWater.txt.in',...
            'Scorer_RBE_libamtrack','TOPAS_scorer_doseRBE_libamtrack.txt.in',...
            'Scorer_RBE_LEM1','TOPAS_scorer_doseRBE_LEM1.txt.in',...
            'Scorer_RBE_WED','TOPAS_scorer_doseRBE_Wedenberg.txt.in',...
            'Scorer_RBE_MCN','TOPAS_scorer_doseRBE_McNamara.txt.in', ...
            ... %PhaseSpace Source
            'phaseSpaceSourcePhotons' ,'VarianClinaciX_6MV_20x20_aboveMLC_w2' );
       

    end

    properties (SetAccess = protected, GetAccess = private)
        topasFolder;
        MCparam; %Struct with parameters of last simulation to be saved to file

        ctR;    %resmpaled CT
    end

    methods
        function obj = matRad_TopasMCEngine(pln)

            if nargin < 1
                pln = [];
            end

            % call superclass constructor
            obj = obj@DoseEngines.matRad_MonteCarloEngineAbstract(pln);
        end

        function setDefaults(this)
            this.setDefaults@DoseEngines.matRad_MonteCarloEngineAbstract();
            matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class

            % Default execution paths are set here
            this.topasFolder = [matRad_cfg.matRadSrcRoot filesep 'doseCalc' filesep 'topas' filesep];
            this.workingDir = [matRad_cfg.primaryUserFolder filesep 'TOPAS' filesep];

            if ~exist(this.workingDir,'dir')
                mkdir(this.workingDir);
                matRad_cfg.dispInfo('Created TOPAS working directory in userfolder %s\n',this.workingDir);
            end
            

            %Let's set some default commands taken from topas installation
            %instructions for mac & debain/ubuntu
            if ispc %We assume topas is installed in wsl (since no windows version)
                this.topasExecCommand = 'wsl export TOPAS_G4_DATA_DIR=~/G4Data; ~/topas/bin/topas';
            elseif ismac
                this.topasExecCommand = 'export TOPAS_G4_DATA_DIR=/Applications/G4Data; export QT_QPA_PLATFORM_PLUGIN_PATH=/Applications/topas/Frameworks; /Applications/topas/bin/topas';
            elseif isunix
                this.topasExecCommand = 'export TOPAS_G4_DATA_DIR=~/G4Data; ~/topas/bin/topas';
            else
                this.topasExecCommand = '';
            end
        end

        function writeAllFiles(obj,ct,cst,stf,machine,w)
            % constructor to write all TOPAS fils for local or external simulation
            %
            % call
            %   topasConfig.writeAllFiles(ct,pln,stf,machine,w)
            %
            % input
            %   ct:             Path to folder where TOPAS files are in (as string)
            %   pln:            matRad plan struct
            %   stf:            matRad steering struct
            %   machine:        machine to be used for calculation
            %   w:              (optional) weights in case of calcDoseDirect

            matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class

            % prepare biological parameters
            if ~isempty(obj.prescribedDose)
                obj.bioParameters.PrescribedDose = obj.prescribedDose;
            end
            if isempty(obj.radiationMode)
                obj.radiationMode = machine.meta.radiationMode;
            end

            % Set correct RBE scorer parameters
            if obj.scorer.RBE
                obj.scorer.doseToMedium = true;
                if any(cellfun(@(teststr) ~isempty(strfind(lower(teststr),'default')), obj.scorer.RBE_model))
                    switch obj.radiationMode
                        case 'protons'
                            obj.scorer.RBE_model = obj.scorer.defaultModelProtons;
                        case {'carbon','helium'}
                            obj.scorer.RBE_model = obj.scorer.defaultModelCarbon;
                        otherwise
                            matRad_cfg.dispError(['No RBE model implemented for ',obj.radiationMode]);
                    end
                end

                % Get alpha beta parameters from bioParam struct
                for i = 1:length(obj.bioParameters.AvailableAlphaXBetaX)
                    if ~isempty(strfind(lower(obj.bioParameters.AvailableAlphaXBetaX{i,2}),'default'))
                        break
                    end
                end
                obj.bioParameters.AlphaX = obj.bioParameters.AvailableAlphaXBetaX{5,1}(1);
                obj.bioParameters.BetaX = obj.bioParameters.AvailableAlphaXBetaX{5,1}(2);

            end
            if obj.scorer.LET
                obj.scorer.doseToMedium = true;
            end

            % create TOPAS working directory if not set
            if ~exist(obj.workingDir,'dir')
                mkdir(obj.workingDir);
                matRad_cfg.dispInfo('Created TOPAS working directory %s\n',obj.workingDir);
            end

            % Write CT, patient parameters and Schneider converter
            matRad_cfg.dispInfo('Writing parameter files to %s\n',obj.workingDir);
            obj.writePatient(ct);

            % Generate uniform weights in case of dij calculation (for later optimization)
            if ~exist('w','var')
                numBixels = sum([stf(:).totalNumOfBixels]);
                w = ones(numBixels,1);
            end

            % Set MCparam structure with important simulation parameters that is needed for later readOut and
            % postprocessing
            obj.MCparam = struct();
            obj.MCparam.tallies = {};
            obj.MCparam.nbRuns = obj.numOfRuns;
            obj.MCparam.simLabel = obj.label;
            obj.MCparam.scoreReportQuantity = obj.scorer.reportQuantity;
            obj.MCparam.workingDir = obj.workingDir;
            obj.MCparam.weights = w;
            obj.MCparam.ctGrid = ct.ctGrid;
            if isfield(ct,'originalGrid')
                obj.MCparam.originalGrid = ct.originalGrid;
            end
            obj.MCparam.cubeDim = ct.cubeDim;
            obj.MCparam.ctResolution = ct.resolution;
            obj.MCparam.numOfCtScen = ct.numOfCtScen;
            % Save used RBE models
            if obj.scorer.RBE
                obj.MCparam.RBE_models = obj.scorer.RBE_model;
                [obj.MCparam.ax,obj.MCparam.bx] = matRad_getPhotonLQMParameters(cst,prod(ct.cubeDim),obj.MCparam.numOfCtScen);
                obj.MCparam.abx(obj.MCparam.bx>0) = obj.MCparam.ax(obj.MCparam.bx>0)./obj.MCparam.bx(obj.MCparam.bx>0);
            end

            % fill in bixels, rays and beams in case of dij calculation or external calculation
            if obj.scorer.calcDij
                counter = 1;
                for f = 1:length(stf)
                    for r = 1:stf(f).numOfRays
                        for b = 1:stf(f).numOfBixelsPerRay(r)
                            obj.MCparam.bixelNum(counter) = b;
                            obj.MCparam.rayNum(counter)   = r;
                            obj.MCparam.beamNum(counter)  = f;
                            counter = counter + 1;
                        end
                    end
                end
            else
                % In case of calcDoseDirect, you only need beamNum
                obj.MCparam.bixelNum = 1;
                obj.MCparam.rayNum   = 1;
                obj.MCparam.beamNum  = 1:length(stf);
            end
            obj.MCparam.numOfRaysPerBeam   = [stf(:).numOfRays];

            % Generate baseData using the MCemittanceBaseData constructor
            % Write TOPAS beam properties
            if ~strcmp(machine.meta.radiationMode,'photons')
                topasBaseData = matRad_MCemittanceBaseData(machine,stf);
            else
                topasBaseData = [];
            end
            obj.writeStfFields(ct,stf,w,topasBaseData);

            % Save simulation parameters to folder
            obj.writeMCparam();

            % Console message
            matRad_cfg.dispInfo('Successfully written TOPAS setup files!\n')
        end

        function dij = readFiles(obj,folder)
            % function to read out TOPAS data
            %
            % call
            %   topasCube = topasConfig.readFiles(folder,dij)
            %   topasCube = obj.readFiles(folder,dij)
            %
            % input
            %   folder:         Path to folder where TOPAS files are in (as string)
            %   dij:            dij struct (this part needs update)
            %
            % output
            %   topasCube:      struct with all read out subfields


            % Load in saved MC parameters
            if isfile([folder filesep 'MCparam.mat'])
                obj.MCparam = load([folder filesep 'MCparam.mat'],'MCparam');
                obj.MCparam = obj.MCparam.MCparam;
            end

            % Read out all TOPAS fields
            topasCubes = obj.readTopasCubes(folder);

            % Set 0 for empty or NaN fields
            topasCubes = obj.markFieldsAsEmpty(topasCubes);

            %% Fill Dij analogously to matRad
            % Prepare empty Dij with empty sparse matrices for fields in topasCubes
            dij = obj.prepareDij(topasCubes);

            % Fill empty Dij with fields from topasCubes
            dij = obj.fillDij(topasCubes,dij);

        end

        function resultGUI = getResultGUI(obj,dij)
            if obj.scorer.calcDij
                resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels,1),dij,1);
            else
                resultGUI = matRad_calcCubes(ones(dij.numOfBeams,1),dij,1);
            end

            % Export RBE model if filled
            if isfield(resultGUI,'RBE_model') && ~isempty(resultGUI.RBE_model)
                resultGUI.RBE_model = dij.RBE_model;
            end

            % Export histories to resultGUI
            if isfield(dij,'nbHistoriesTotal')
                resultGUI.nbHistoriesTotal = dij.nbHistoriesTotal;
                resultGUI.nbParticlesTotal = dij.nbParticlesTotal;
            end

            % stuff for 4D
            %             if this.multScen.totNumScen ~= 1
            %                 resultGUI.accPhysicalDose = zeros(size(resultGUI.phaseDose{1}));
            %                 for i = 1:this.multScen.totNumScen
            %                     resultGUI.accPhysicalDose = resultGUI.accPhysicalDose + resultGUI.phaseDose{i};
            %                 end
            %             end
        end

        function dij = readExternal(obj,folder)
            % function to read out complete TOPAS simulation from single folder
            %
            % call
            %   topasCube = topasConfig.readExternal(folder)
            %   topasCube = obj.readExternal(folder)
            %
            % input
            %   folder:         Path to folder where TOPAS files are in (as string)
            %
            % output
            %   topasCube:      struct with all read out subfields
            %
            % EXAMPLE calls:
            %   topasCube = topasConfig.readExternal('pathToFolder')

            % read in TOPAS files in dij
            dij = obj.readFiles(folder);

            % Postprocessing
            %resultGUI = obj.getResultGUI(dij);

        end
    end

    methods (Access = protected)
        function dij = calcDose(this,ct,cst,stf)
            % Instance of MatRad_Config class
            matRad_cfg = MatRad_Config.instance();

            % Set parameters for full Dij calculation
            if ~this.calcDoseDirect
                this.scorer.calcDij = true;
                this.numOfRuns = 1;
            end
            % set nested folder structure if external calculation is turned on (this will put new simulations in subfolders)
            if strcmp(this.externalCalculation,'write')
                this.workingDir = [matRad_cfg.primaryUserFolder filesep 'TOPAS' filesep];
                this.workingDir = [this.workingDir stf(1).radiationMode,'_',stf(1).machine,'_',datestr(now, 'dd-mm-yy')];
            elseif isfolder(this.externalCalculation)
                dij = this.readExternal(this.externalCalculation);
                return;
            else
            end
                

            %% Initialize dose grid and dij

            % load calcDoseInit as usual
            dij = this.initDoseCalc(ct,cst,stf);

            %% sending data to topas
            if isfield(this.machine.meta,'SCD')
                for i = 1:size(stf,2)
                    stf(i).SCD = this.machine.meta.SCD;
                end
            end

            % Collect given weights
            if this.calcDoseDirect
                %     w = zeros(sum([stf(:).totalNumOfBixels]),ctR.numOfCtScen);
                w = zeros(sum([stf(:).totalNumOfBixels]),1);
                counter = 1;
                for i = 1:length(stf)
                    for j = 1:stf(i).numOfRays
                        rayBix = stf(i).numOfBixelsPerRay(j);
                        if isfield(stf(1).ray, 'shapes')
                            w(counter:counter+rayBix-1)  = [stf(i).ray.shapes.weight];
                        else
                            w(counter:counter+rayBix-1,:) = stf(i).ray(j).weight;
                        end
                        counter = counter + rayBix;
                    end
                end
            end
            
            for i = 1:numel(stf)
                if strcmp(stf(i).radiationMode,'photons')
                    stf(i).ray.energy = stf(i).ray.energy.*ones(size(w));
                end
            end            

            % Get photon parameters for RBExD calculation
            if this.calcBioDose
                this.scorer.RBE = true;
                [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,1,VdoseGrid);
                dij.abx(dij.bx>0) = dij.ax(dij.bx>0)./dij.bx(dij.bx>0);
            end

            % save current directory to revert back to later
            currDir = cd;

            for shiftScen = 1:this.multScen.totNumShiftScen

                %Find first instance of the shift to select the shift values
                ixShiftScen = find(this.multScen.linearMask(:,2) == shiftScen,1);

                % manipulate isocenter
                for k = 1:numel(stf)
                    stf(k).isoCenter = stf(k).isoCenter + this.multScen.isoShift(ixShiftScen,:);
                end

                % Delete previous topas files so there is no mix-up
                files = dir([this.workingDir,'*']);
                files = {files(~[files.isdir]).name};
                fclose('all');
                for i = 1:length(files)
                    delete([this.workingDir,files{i}])
                end

                % Run simulations for each scenario
                for ctScen = 1:this.multScen.numOfCtScen
                    for rangeShiftScen = 1:this.multScen.totNumRangeScen
                        if this.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)

                            % Save ctScen and rangeShiftScen for file constructor
                            if ct.numOfCtScen > 1
                                this.ctR.currCtScen = ctScen;
                                this.ctR.currRangeShiftScen = rangeShiftScen;
                            end

                            % actually write TOPAS files
                            if this.calcDoseDirect
                                this.writeAllFiles(this.ctR,cst,stf,this.machine,w);
                            else
                                this.writeAllFiles(this.ctR,cst,stf,this.machine);
                            end
                        end
                    end
                end

                % change director back to original directory
                cd(this.workingDir);

                % Skip local calculation and data readout with this parameter. All necessary parameters to read the data back in
                % later are stored in the MCparam file that is stored in the folder. The folder is generated in the working
                % directory and the matRad_plan*.txt file can be manually called with TOPAS.
                if strcmp(this.externalCalculation,'write')
                    matRad_cfg.dispInfo(['TOPAS simulation skipped for external calculation\nFiles have been written to: "',replace(this.workingDir,'\','\\'),'"']);
                else
                    for ctScen = 1:ct.numOfCtScen
                        for beamIx = 1:numel(stf)
                            for runIx = 1:this.numOfRuns
                                if ct.numOfCtScen > 1
                                    fname = sprintf('%s_field%d_ct%d_run%d',this.label,beamIx,ctScen,runIx);
                                else
                                    fname = sprintf('%s_field%d_run%d',this.label,beamIx,runIx);
                                end

                                if isprop(this,'verbosity') && strcmp(this.verbosity,'full')
                                    topasCall = sprintf('%s %s.txt',this.topasExecCommand,fname);
                                else
                                    topasCall = sprintf('%s %s.txt > %s.out > %s.log',this.topasExecCommand,fname,fname,fname);
                                end

                                % initiate parallel runs and delete previous files
                                if this.parallelRuns
                                    finishedFiles{runIx} = sprintf('%s.finished',fname);
                                    topasCall = [topasCall '; touch ' finishedFiles{runIx} ' &'];
                                end

                                % Actual simulation happening here
                                matRad_cfg.dispInfo('Calling TOPAS: %s\n',topasCall);
                                [status,cmdout] = system(topasCall,'-echo');

                                % Process TOPAS output and potential errors
                                cout = splitlines(string(cmdout));
                                if status == 0
                                    matRad_cfg.dispInfo('TOPAS simulation completed succesfully\n');
                                else
                                    if status == 139
                                        matRad_cfg.dispError('TOPAS segmentation fault: might be caused from an outdated TOPAS version or Linux distribution');
                                    else
                                        matRad_cfg.dispError('TOPAS simulation exited with error code %d\n "%s"',status,cout(2:end-1));
                                    end
                                end
                            end

                            % wait for parallel runs to finish and process
                            if this.parallelRuns
                                runsFinished = false;
                                pause('on');
                                while ~runsFinished
                                    pause(1);
                                    fin = cellfun(@(f) exist(f,'file'),finishedFiles);
                                    runsFinished = all(fin);
                                end
                                % Delete marker files
                                delete(finishedFiles{:});
                            end
                        end
                    end
                end

                % revert back to original directory
                cd(currDir);

                % manipulate isocenter back
                for k = 1:length(stf)
                    stf(k).isoCenter = stf(k).isoCenter - this.multScen.isoShift(ixShiftScen,:);
                end

            end

            %% Simulation(s) finished - read out volume scorers from topas simulation
            % Skip readout if external files were generated
            if strcmp(this.externalCalculation,'off')
                dij = this.readFiles(this.workingDir);

                % Order fields for easier comparison between different dijs
                dij = orderfields(dij);
            else
                dij.beamNum = 1;
                dij.bixelNum = 1;
                dij.ctGrid = this.ctR.ctGrid;
                dij.doseGrid = this.doseGrid;                
                dij.numOfBeams = 1;
                dij.numOfRaysPerBeam = 1;
                dij.numOfScenarios = this.multScen.totNumScen;
                for i = 1:this.multScen.numOfCtScen
                    for j = 1:this.multScen.totNumShiftScen
                        for k = 1:this.multScen.totNumRangeScen
                            if this.multScen.scenMask(i,j,k)
                                %TODO: loop over all expected output quantities
                                dij.physicalDose{i,j,k} = zeros(dij.ctGrid.numOfVoxels,1);
                                dij.physicalDose_std{i,j,k} = zeros(dij.ctGrid.numOfVoxels,1);
                            end

                        end
                    end
                end
                dij.rayNum = 1;
                dij.totalNumOfBixels = 1;
                dij.totalNumOfRays = 1;
                dij.meta.TOPASworkingDir = this.workingDir;
            end

            this.finalizeDose();

        end

                
        function dij = initDoseCalc(this,ct,cst,stf)
            dij = this.initDoseCalc@DoseEngines.matRad_MonteCarloEngineAbstract(ct,cst,stf);
            matRad_cfg = MatRad_Config.instance();

            % % for TOPAS we explicitly downsample the ct to the dose grid (might not be necessary in future versions with separated grids)
            % Check if CT has already been resampled
            matRad_cfg.dispInfo('Resampling cst... ');
            if ~isfield(ct,'resampled')
                % Allpcate resampled cubes
                cubeHUresampled = cell(1,ct.numOfCtScen);
                cubeResampled = cell(1,ct.numOfCtScen);
            
                % Perform resampling to dose grid
                for s = 1:ct.numOfCtScen
                    cubeHUresampled{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cubeHU{s}, ...
                        dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
                    cubeResampled{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cube{s}, ...
                        dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
                end
            
                % Allocate temporary resampled CT
                this.ctR = ct;
                this.ctR.cube = cell(1);
                this.ctR.cubeHU = cell(1);
            
                % Set CT resolution to doseGrid resolution
                this.ctR.resolution = dij.doseGrid.resolution;
                this.ctR.cubeDim = dij.doseGrid.dimensions;
                this.ctR.x = dij.doseGrid.x;
                this.ctR.y = dij.doseGrid.y;
                this.ctR.z = dij.doseGrid.z;
            
                % Write resampled cubes
                this.ctR.cubeHU = cubeHUresampled;
                this.ctR.cube = cubeResampled;
            
                % Set flag for complete resampling
                this.ctR.resampled = 1;
                this.ctR.ctGrid = dij.doseGrid;
            
                % Save original grid
                this.ctR.originalGrid = dij.ctGrid;
                matRad_cfg.dispInfo('done!\n');
            else
                this.ctR = ct;
                matRad_cfg.dispInfo('already resampled. Skipping! \n');
            end

            % overwrite CT grid in dij in case of modulation.
            if isfield(this.ctR,'ctGrid')
                dij.ctGrid = this.ctR.ctGrid;
            end
        end
    end
    methods (Access = private)
        function topasCubes = markFieldsAsEmpty(obj,topasCubes)

            matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class

            % Check if all fields in topasCubes are filled or overwrite with 0 if not.
            fields = fieldnames(topasCubes);
            for field = 1:length(fields)
                if all(isnan(topasCubes.(fields{field}){1}(:)) | topasCubes.(fields{field}){1}(:)==0)
                    matRad_cfg.dispWarning(['Field ' fields{field} ' in topasCubes resulted in all zeros and NaN.'])
                    topasCubes.(fields{field}) = 0;
                end
            end

        end

        function topasCube = readTopasCubes(obj,folder)
            % function to read out TOPAS data
            %
            % call
            %   topasCube = topasConfig.readTopasCubes(folder,dij)
            %   topasCube = obj.readTopasCubes(folder,dij)
            %
            % input
            %   folder:         Path to folder where TOPAS files are in (as string)
            %   dij:            dij struct (this part needs update)
            %
            % output
            %   topasCube:      struct with all read out subfields

            matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class

            %          obj.MCparam.scoreReportQuantity = 'Sum';
            % Process reportQuantities (for example 'Sum' or 'Standard_Deviation'read
            if iscell(obj.MCparam.scoreReportQuantity)
                obj.MCparam.numOfReportQuantities = length(obj.MCparam.scoreReportQuantity);
            else
                obj.MCparam.numOfReportQuantities = 1;
                obj.MCparam.scoreReportQuantity = {obj.MCparam.scoreReportQuantity};
            end

            % Normalize with histories and particles/weight
            correctionFactor = obj.numParticlesPerHistory / double(obj.MCparam.nbHistoriesTotal);

            % Get all saved quantities
            % Make sure that the filename always ends on 'run1_tally'
            switch obj.MCparam.outputType
                case 'csv'
                    searchstr = 'score_matRad_plan_field1_run1_*.csv';
                    files = dir([folder filesep searchstr]);
                    %obj.MCparam.tallies = cellfun(@(x) extractBetween(x,'run1_','.csv') ,{files(:).name}); %Not Octave compatible
                    nameBegin = strfind(searchstr,'*');
                    obj.MCparam.tallies = cellfun(@(s) s(nameBegin:end-4),{files(:).name},'UniformOutput',false);
                case 'binary'
                    searchstr = 'score_matRad_plan_field1_run1_*.bin';
                    files = dir([folder filesep searchstr]);
                    %obj.MCparam.tallies = cellfun(@(x) extractBetween(x,'run1_','.bin') ,{files(:).name}); %Not Octave compatible
                    nameBegin = strfind(searchstr,'*');
                    obj.MCparam.tallies = cellfun(@(s) s(nameBegin:end-4),{files(:).name},'UniformOutput',false);
            end

            obj.MCparam.tallies = unique(obj.MCparam.tallies);
            talliesCut = replace(obj.MCparam.tallies,'-','_');

            % Load data for each tally individually
            for t = 1:length(obj.MCparam.tallies)
                tnameFile = obj.MCparam.tallies{t};
                tname = talliesCut{t};
                % Loop over all beams/fields and ctScenarios
                for f = 1:obj.MCparam.nbFields
                    for ctScen = 1:obj.MCparam.numOfCtScen

                        % Loop over all batches/runs
                        for k = 1:obj.MCparam.nbRuns
                            % Get file name of current field, run and tally (and ct, if applicable)
                            if obj.MCparam.numOfCtScen > 1
                                genFileName = sprintf('score_%s_field%d_ct%d_run%d_%s',obj.MCparam.simLabel,f,ctScen,k,tnameFile);
                            else
                                genFileName = sprintf('score_%s_field%d_run%d_%s',obj.MCparam.simLabel,f,k,tnameFile);
                            end


                            switch obj.MCparam.outputType
                                case 'csv'
                                    % Generate csv file path to load
                                    genFullFile = fullfile(folder,[genFileName '.csv']);
                                case 'binary'
                                    % Generate bin file path to load
                                    genFullFile = fullfile(folder,[genFileName '.bin']);
                                otherwise
                                    matRad_cfg.dispError('Not implemented!');
                            end

                            % Read data from scored TOPAS files
                            dataRead = obj.readBinCsvData(genFullFile);

                            for i = 1:numel(dataRead)
                                data.(obj.MCparam.scoreReportQuantity{i}){k} = dataRead{i};
                            end

                            % for example the standard deviation is not calculated for alpha/beta so a loop through all
                            % reportQuantities does not work here
                            currNumOfQuantities = numel(dataRead);
                        end

                        % Set dimensions of output cube
                        cubeDim = size(dataRead{1});

                        % add STD quadratically
                        for i = 1:currNumOfQuantities
                            if ~isempty(strfind(lower(obj.MCparam.scoreReportQuantity{i}),'standard_deviation'))
                                topasSum.(obj.MCparam.scoreReportQuantity{i}) = sqrt(double(obj.MCparam.nbHistoriesTotal)) * sqrt(sum(cat(4,data.(obj.MCparam.scoreReportQuantity{i}){:}).^2,4));
                            else
                                topasSum.(obj.MCparam.scoreReportQuantity{i}) = sum(cat(4,data.(obj.MCparam.scoreReportQuantity{i}){:}),4);
                            end
                        end

                        if ~isempty(strfind(lower(tnameFile),'dose'))
                            if obj.MCparam.nbRuns > 1
                                % Calculate Standard Deviation from batches
                                topasMeanDiff = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
                                for k = 1:obj.MCparam.nbRuns
                                    topasMeanDiff = topasMeanDiff + (data.Sum{k} - topasSum.Sum / obj.MCparam.nbRuns).^2;
                                end
                                % variance of the mean
                                topasVarMean = topasMeanDiff./(obj.MCparam.nbRuns - 1)./obj.MCparam.nbRuns;
                                % std of the MEAN!
                                topasStdMean = sqrt(topasVarMean);
                                % std of the SUM
                                topasStdSum = topasStdMean * correctionFactor * obj.MCparam.nbRuns;

                                % Save std to topasCube
                                topasCube.([tname '_batchStd_beam' num2str(f)]){ctScen} = topasStdSum;
                            end

                            for i = 1:currNumOfQuantities
                                topasSum.(obj.MCparam.scoreReportQuantity{i}) = correctionFactor .* topasSum.(obj.MCparam.scoreReportQuantity{i});
                            end

                        elseif any(cellfun(@(teststr) ~isempty(strfind(tname,teststr)), {'alpha','beta','RBE','LET'}))
                            for i = 1:currNumOfQuantities
                                topasSum.(obj.MCparam.scoreReportQuantity{i}) = topasSum.(obj.MCparam.scoreReportQuantity{i}) ./ obj.MCparam.nbRuns;
                            end
                        end

                        % Tally per field
                        if isfield(topasSum,'Sum')
                            topasCube.([tname '_beam' num2str(f)]){ctScen} = topasSum.Sum;
                        end
                        if isfield(topasSum,'Standard_Deviation')
                            topasCube.([tname '_std_beam' num2str(f)]){ctScen} = topasSum.Standard_Deviation;
                        end
                    end
                end
            end
        end

        function dataOut = readBinCsvData(~,genFullFile)

            matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class

            if ~isempty(strfind(lower(genFullFile),'.csv'))
                % Read csv header file to get cubeDim and number of scorers automatically
                fid = fopen(genFullFile);
                header = textscan(fid,'%[^,],%[^,],%[^,]',1);
                fclose(fid);

                % Split header in rows
                header = strsplit(strrep(header{1}{1},' ',''),'#');
            elseif ~isempty(strfind(lower(genFullFile),'.bin'))
                % Isolate filename without ending
                [folder, filename] = fileparts(genFullFile);
                strippedFileName = [folder filesep filename];

                % Read binheader file to get cubeDim and number of scorers automatically
                fID = fopen([strippedFileName '.binheader']);
                header = textscan(fID,'%c');
                fclose(fID);

                % Split header in rows
                header = strsplit(header{1}','#')';
            else
                % Error if neither csv nor bin
                matRad_cfg.dispError('Not implemented!');
            end

            % Find rows where number of bins are stored
            xLine = find(cellfun(@(x) ~isempty(x), strfind(header,'Xin')));
            cubeDim(2) = str2double(header{xLine}(4:strfind(header{xLine},'binsof')-1));
            cubeDim(1) = str2double(header{xLine+1}(4:strfind(header{xLine+1},'binsof')-1));
            cubeDim(3) = str2double(header{xLine+2}(4:strfind(header{xLine+2},'binsof')-1));

            if ~isempty(strfind(lower(genFullFile),'.csv'))
                % Read out bin data
                dataOut = matRad_readCsvData(genFullFile,cubeDim);
            elseif ~isempty(strfind(lower(genFullFile),'.bin'))
                % Read out bin data
                dataOut = matRad_readBinData(genFullFile,cubeDim);
            end

        end

        function dij = prepareDij(obj,topasCubes)

            % Load ctScen variable
            numOfScenarios = obj.MCparam.numOfCtScen;

            % Set flag for RBE and LET
            if any(cellfun(@(teststr) ~isempty(strfind(lower(teststr),'alpha')), fieldnames(topasCubes)))
                obj.scorer.RBE = true;
            end
            if any(cellfun(@(teststr) ~isempty(strfind(teststr,'LET')), fieldnames(topasCubes)))
                obj.scorer.LET = true;
            end

            % Create empty dij
            dij.numOfScenarios = numOfScenarios;
            dij.numOfBeams = max(obj.MCparam.beamNum);
            dij.numOfRaysPerBeam = obj.MCparam.numOfRaysPerBeam;
            dij.totalNumOfRays = sum(dij.numOfRaysPerBeam);
            dij.totalNumOfBixels = length(obj.MCparam.bixelNum);
            dij.bixelNum = obj.MCparam.bixelNum';
            dij.rayNum = obj.MCparam.rayNum';
            dij.beamNum = obj.MCparam.beamNum';

            % Write dij grids
            dij.doseGrid = obj.MCparam.ctGrid;
            dij.ctGrid = obj.MCparam.originalGrid;

            % Save RBE models in dij for postprocessing in calcCubes
            if obj.scorer.RBE
                dij.RBE_models = obj.MCparam.RBE_models;
                dij.ax = obj.MCparam.ax;
                dij.bx = obj.MCparam.bx;
                dij.abx = obj.MCparam.abx;
            end

            % Get basic tallies from topasCubes for sparse matrix allocation
            beamNames = strsplit(sprintf('_beam%i,',1:dij.numOfBeams),',');
            if obj.scorer.calcDij
                rayNames = strsplit(sprintf('_ray%i,',unique(dij.rayNum)),',');
                bixelNames = strsplit(sprintf('_bixel%i,',unique(dij.bixelNum)),',');
                topasCubesTallies = unique(erase(fieldnames(topasCubes),rayNames));
                topasCubesTallies = unique(erase(topasCubesTallies,bixelNames));
                topasCubesTallies = unique(erase(topasCubesTallies,beamNames));
            else
                topasCubesTallies = unique(erase(fieldnames(topasCubes),beamNames));
            end

            % Get default tallies from dose
            dijTallies = topasCubesTallies(cellfun(@(teststr) ~isempty(strfind(lower(teststr),'dose')), topasCubesTallies));

            % Handle LET tally
            if obj.scorer.LET
                dijTallies{end+1} = 'mLETDose';
                %                 dijTallies{end+1} = 'LET';
            end

            % Get unique tallies for RBE models
            if obj.scorer.RBE
                for r = 1:length(obj.MCparam.RBE_models)
                    dijTallies{end+1} = ['mAlphaDose_' obj.MCparam.RBE_models{r}];
                    dijTallies{end+1} = ['mSqrtBetaDose_' obj.MCparam.RBE_models{r}];
                    %                     dijTallies{end+1} = 'alpha';
                    %                     dijTallies{end+1} = 'beta';
                end
            end

            % Create empty sparse matrices
            % Note that for MonteCarlo, there are no individual bixels, but only 2 beams
            for t = 1:length(dijTallies)
                for ctScen = 1:dij.numOfScenarios
                    if obj.scorer.calcDij
                        dij.(dijTallies{t}){ctScen,1} = spalloc(obj.MCparam.ctGrid.numOfVoxels,dij.totalNumOfBixels,1);
                    else
                        dij.(dijTallies{t}){ctScen,1} = spalloc(obj.MCparam.ctGrid.numOfVoxels,dij.numOfBeams,1);
                    end
                end
            end

        end

        function dij = fillDij(obj,topasCubes,dij)
            %TODO: Insert documentation
            matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class

            % Load weights from parameter variable
            w = obj.MCparam.weights;

            % Get basic tallies from topasCubes for sparse matrix allocation
            beamNames = strsplit(sprintf('_beam%i,',1:dij.numOfBeams),',');
            if obj.scorer.calcDij
                rayNames = strsplit(sprintf('_ray%i,',unique(dij.rayNum)),',');
                bixelNames = strsplit(sprintf('_bixel%i,',unique(dij.bixelNum)),',');
                topasCubesTallies = unique(erase(fieldnames(topasCubes),rayNames));
                topasCubesTallies = unique(erase(topasCubesTallies,bixelNames));
                topasCubesTallies = unique(erase(topasCubesTallies,beamNames));
            else
                topasCubesTallies = unique(erase(fieldnames(topasCubes),beamNames));
            end

            % Allocate possible scored quantities
            processedQuantities = {'','_std','_batchStd'};
            topasCubesTallies = unique(erase(topasCubesTallies,processedQuantities(2:end)));
            

            % Loop through 4D scenarios
            for ctScen = 1:dij.numOfScenarios

                % Process physicalDose
                % this is done separately since it's needed for processing the other dose fields
                if obj.scorer.calcDij
                    for d = 1:dij.totalNumOfBixels
                        physDoseFields = strfind(lower(topasCubesTallies),'physicaldose');
                        physDoseFields = not(cellfun('isempty',physDoseFields));
                        for j = find(physDoseFields)'
                            % loop through possible quantities
                            for p = 1:length(processedQuantities)
                                % Check if current quantity is available and write to dij
                                if isfield(topasCubes,[topasCubesTallies{j} '_ray' num2str(dij.rayNum(d)) '_bixel' num2str(dij.bixelNum(d)) processedQuantities{p} '_beam' num2str(dij.beamNum(d))]) ...
                                        && iscell(topasCubes.([topasCubesTallies{j} '_ray' num2str(dij.rayNum(d)) '_bixel' num2str(dij.bixelNum(d)) processedQuantities{p} '_beam' num2str(dij.beamNum(d))]))
                                    dij.([topasCubesTallies{j} processedQuantities{p}]){ctScen,1}(:,d) = sum(w)*reshape(topasCubes.([topasCubesTallies{j} '_ray' num2str(dij.rayNum(d)) '_bixel' num2str(dij.bixelNum(d)) processedQuantities{p} '_beam' num2str(dij.beamNum(d))]){ctScen},[],1);
                                end
                            end
                        end
                    end
                else
                    for d = 1:dij.numOfBeams
                        physDoseFields = strfind(lower(topasCubesTallies),'physicaldose');
                        physDoseFields = not(cellfun('isempty',physDoseFields));
                        for j = find(physDoseFields)'
                            for p = 1:length(processedQuantities)
                                % Check if current quantity is available and write to dij
                                if isfield(topasCubes,[topasCubesTallies{j} processedQuantities{p} '_beam' num2str(d)]) && iscell(topasCubes.([topasCubesTallies{j} processedQuantities{p} '_beam' num2str(d)]))
                                    dij.([topasCubesTallies{j} processedQuantities{p}]){ctScen}(:,d) = sum(w)*reshape(topasCubes.([topasCubesTallies{j} processedQuantities{p} '_beam',num2str(d)]){ctScen},[],1);
                                end
                            end
                        end
                    end
                end

                % Remove processed physDoseFields from total tallies
                topasCubesTallies = topasCubesTallies(~physDoseFields);

                % Process other fields
                if obj.scorer.calcDij
                    for d = 1:dij.totalNumOfBixels
                        for j = 1:numel(topasCubesTallies)
                            % Handle dose to water
                            if ~isempty(strfind(lower(topasCubesTallies{j}),'dose'))
                                % loop through possible quantities
                                for p = 1:length(processedQuantities)
                                    % Check if current quantity is available and write to dij
                                    if isfield(topasCubes,[topasCubesTallies{j} '_ray' num2str(dij.rayNum(d)) '_bixel' num2str(dij.bixelNum(d)) processedQuantities{p} '_beam' num2str(dij.beamNum(d))]) ...
                                            && iscell(topasCubes.([topasCubesTallies{j} '_ray' num2str(dij.rayNum(d)) '_bixel' num2str(dij.bixelNum(d)) processedQuantities{p} '_beam' num2str(dij.beamNum(d))]))
                                        dij.([topasCubesTallies{j} processedQuantities{p}]){ctScen,1}(:,d) = sum(w)*reshape(topasCubes.([topasCubesTallies{j} '_ray' num2str(dij.rayNum(d)) '_bixel' num2str(dij.bixelNum(d)) processedQuantities{p} '_beam' num2str(dij.beamNum(d))]){ctScen},[],1);
                                    end
                                end
                                % Handle RBE-related quantities (not multiplied by sum(w)!)
                            elseif ~isempty(strfind(lower(topasCubesTallies{j}),'alpha'))
                                modelName = strsplit(topasCubesTallies{j},'_');
                                modelName = modelName{end};
                                if isfield(topasCubes,[topasCubesTallies{j} '_ray' num2str(dij.rayNum(d)) '_bixel' num2str(dij.bixelNum(d)) '_beam' num2str(dij.beamNum(d))]) ...
                                        && iscell(topasCubes.([topasCubesTallies{j} '_ray' num2str(dij.rayNum(d)) '_bixel' num2str(dij.bixelNum(d)) '_beam' num2str(dij.beamNum(d))]))
                                    dij.(['mAlphaDose_' modelName]){ctScen,1}(:,d) = reshape(topasCubes.([topasCubesTallies{j} '_ray' num2str(dij.rayNum(d)) '_bixel' num2str(dij.bixelNum(d)) '_beam' num2str(dij.beamNum(d))]){ctScen},[],1) .* dij.physicalDose{ctScen,1}(:,d);
                                end
                            elseif ~isempty(strfind(lower(topasCubesTallies{j}),'beta'))
                                modelName = strsplit(topasCubesTallies{j},'_');
                                modelName = modelName{end};
                                if isfield(topasCubes,[topasCubesTallies{j} '_ray' num2str(dij.rayNum(d)) '_bixel' num2str(dij.bixelNum(d)) '_beam' num2str(dij.beamNum(d))]) ...
                                        && iscell(topasCubes.([topasCubesTallies{j} '_ray' num2str(dij.rayNum(d)) '_bixel' num2str(dij.bixelNum(d)) '_beam' num2str(dij.beamNum(d))]))
                                    dij.(['mSqrtBetaDose_' modelName]){ctScen,1}(:,d) = sqrt(reshape(topasCubes.([topasCubesTallies{j} '_ray' num2str(dij.rayNum(d)) '_bixel' num2str(dij.bixelNum(d)) '_beam' num2str(dij.beamNum(d))]){ctScen},[],1)) .* dij.physicalDose{ctScen,1}(:,d);
                                end
                            elseif ~isempty(strfind(topasCubesTallies{j},'LET'))
                                if isfield(topasCubes,[topasCubesTallies{j} '_ray' num2str(dij.rayNum(d)) '_bixel' num2str(dij.bixelNum(d)) '_beam' num2str(dij.beamNum(d))]) ...
                                        && iscell(topasCubes.([topasCubesTallies{j} '_ray' num2str(dij.rayNum(d)) '_bixel' num2str(dij.bixelNum(d)) '_beam' num2str(dij.beamNum(d))]))
                                    dij.mLETDose{ctScen,1}(:,d) = reshape(topasCubes.([topasCubesTallies{j} '_ray' num2str(dij.rayNum(d)) '_bixel' num2str(dij.bixelNum(d)) '_beam' num2str(dij.beamNum(d))]){ctScen},[],1) .* dij.physicalDose{ctScen,1}(:,d);
                                end
                            else
                                matRad_cfg.dispError('Postprocessing error: Tallies handles incorrectly')
                            end
                        end
                    end
                else
                    for d = 1:dij.numOfBeams
                        for j = 1:numel(topasCubesTallies)
                            % Handle dose to medium and dose to water
                            if ~isempty(strfind(lower(topasCubesTallies{j}),'dose'))
                                % loop through possible quantities
                                for p = 1:length(processedQuantities)
                                    % Check if current quantity is available and write to dij
                                    if isfield(topasCubes,[topasCubesTallies{j} processedQuantities{p} '_beam' num2str(d)]) && iscell(topasCubes.([topasCubesTallies{j} processedQuantities{p} '_beam' num2str(d)]))
                                        dij.([topasCubesTallies{j} processedQuantities{p}]){ctScen}(:,d) = sum(w)*reshape(topasCubes.([topasCubesTallies{j} processedQuantities{p} '_beam',num2str(d)]){ctScen},[],1);
                                    end
                                end
                                % Handle RBE-related quantities (not multiplied by sum(w)!)
                            elseif ~isempty(strfind(lower(topasCubesTallies{j}),'alpha'))
                                modelName = strsplit(topasCubesTallies{j},'_');
                                modelName = modelName{end};
                                if isfield(topasCubes,[topasCubesTallies{j} '_beam' num2str(d)]) && iscell(topasCubes.([topasCubesTallies{j} '_beam' num2str(d)]))
                                    dij.(['mAlphaDose_' modelName]){ctScen}(:,d)        = reshape(topasCubes.([topasCubesTallies{j} '_beam',num2str(d)]){ctScen},[],1) .* dij.physicalDose{ctScen}(:,d);
                                end
                            elseif ~isempty(strfind(lower(topasCubesTallies{j}),'beta'))
                                modelName = strsplit(topasCubesTallies{j},'_');
                                modelName = modelName{end};
                                if isfield(topasCubes,[topasCubesTallies{j} '_beam' num2str(d)]) && iscell(topasCubes.([topasCubesTallies{j} '_beam' num2str(d)]))
                                    dij.(['mSqrtBetaDose_' modelName]){ctScen}(:,d)        = sqrt(reshape(topasCubes.([topasCubesTallies{j} '_beam',num2str(d)]){ctScen},[],1)) .* dij.physicalDose{ctScen}(:,d);
                                end
                            elseif ~isempty(strfind(topasCubesTallies{j},'LET'))
                                if isfield(topasCubes,[topasCubesTallies{j} '_beam' num2str(d)]) && iscell(topasCubes.([topasCubesTallies{j} '_beam' num2str(d)]))
                                    dij.mLETDose{ctScen}(:,d)        = reshape(topasCubes.([topasCubesTallies{j} '_beam',num2str(d)]){ctScen},[],1) .* dij.physicalDose{ctScen}(:,d);
                                end
                            else
                                matRad_cfg.dispError('Postprocessing error: Tallies handles incorrectly')
                            end
                        end
                    end
                end
            end

            % Save number of histories and particles to Dij
            dij.nbHistoriesTotal = obj.MCparam.nbHistoriesTotal;
            dij.nbParticlesTotal = obj.MCparam.nbParticlesTotal;

        end

        function writeRunHeader(obj,fID,fieldIx,runIx,ctScen)
            %TODO: Insert documentation
            fprintf(fID,'s:Sim/PlanLabel = "%s"\n',obj.label);
            if exist('ctScen','var')
                fprintf(fID,'s:Sim/ScoreLabel = "score_%s_field%d_ct%d_run%d"\n',obj.label,fieldIx,ctScen,runIx);
            else
                fprintf(fID,'s:Sim/ScoreLabel = "score_%s_field%d_run%d"\n',obj.label,fieldIx,runIx);
            end
            fprintf(fID,'\n');

            logicalString = {'"False"', '"True"'};

            fprintf(fID,'i:Ma/Verbosity = %d\n',obj.verbosity.material);
            fprintf(fID,'i:Ts/TrackingVerbosity = %d\n',obj.verbosity.tracking);
            fprintf(fID,'i:Ts/EventVerbosity = %d\n',obj.verbosity.event);
            fprintf(fID,'i:Ts/RunVerbosity = %d\n',obj.verbosity.run);
            fprintf(fID,'b:Ts/ShowCPUTime = %s\n',logicalString{obj.verbosity.cputime + 1});
            fprintf(fID,'i:Tf/Verbosity = %d\n',obj.verbosity.timefeatures);
            fprintf(fID,'i:Ts/MaxInterruptedHistories = %d\n',obj.verbosity.maxinterruptedhistories);
            fprintf(fID,'i:Ts/NumberOfThreads = %d\n',obj.numThreads);
            fprintf(fID,'i:Ts/MaximumNumberOfDetailedErrorReports = %d\n',obj.verbosity.maxDetailedErrorReports);
            fprintf(fID,'i:Ts/ShowHistoryCountAtInterval = %d\n',10^(floor(log10(1/obj.numOfRuns * obj.numHistoriesDirect))-1));
            fprintf(fID,'\n');


            fprintf(fID,'s:Sim/DoseScorerOutputType = "%s"\n',obj.scorer.outputType);
            if iscell(obj.scorer.reportQuantity)
                fprintf(fID,'sv:Sim/DoseScorerReport = %i ',length(obj.scorer.reportQuantity));
                fprintf(fID,'"%s" ',obj.scorer.reportQuantity{:});
                fprintf(fID,'\n');
            else
                fprintf(fID,'sv:Sim/DoseScorerReport = 1 "%s"\n',obj.scorer.reportQuantity);
            end
            fprintf(fID,'\n');
            fprintf(fID,['i:Ts/Seed = ',num2str(runIx),'\n']);

            %TODO: remove or document 
            %fprintf(fID,'includeFile = %s/TOPAS_Simulation_Setup.txt\n',obj.thisFolder);
            %fprintf(fID,'includeFile = %s/TOPAS_matRad_geometry.txt\n',obj.thisFolder);
            %fprintf(fID,'includeFile = %s/TOPAS_scorer_surfaceIC.txt\n',obj.thisFolder);
        end

        function writeFieldHeader(obj,fID,ctScen)
            %TODO: Insert documentation
            matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class

            if ~strcmp(obj.beamProfile,'phasespace')
                fprintf(fID,'u:Sim/HalfValue = %d\n',0.5);
                fprintf(fID,'u:Sim/SIGMA2FWHM = %d\n',2.354818);
                fprintf(fID,'u:Sim/FWHM2SIGMA = %d\n',0.424661);
                fprintf(fID,'\n');
            else
                fprintf(fID,'u:Sim/HalfValue = %d\n',0.5);
                fprintf(fID,'u:Sim/SIGMA2FWHM = %d\n',1.1905);
                fprintf(fID,'u:Sim/FWHM2SIGMA = %d\n',0.84);
                fprintf(fID,'\n');
            end

            fprintf(fID,'d:Sim/ElectronProductionCut = %f mm\n',obj.electronProductionCut);
            fprintf(fID,'s:Sim/WorldMaterial = "%s"\n',obj.worldMaterial);
            fprintf(fID,'\n');

            % Add ctScen number to filenames
            if exist('ctScen','var')
                paramFile = strsplit(obj.outfilenames.patientParam,'.');
                paramFile = strjoin(paramFile,[num2str(ctScen) '.']);
            else
                paramFile = obj.outfilenames.patientParam;
            end

            fprintf(fID,'includeFile = %s\n',paramFile);
            fprintf(fID,'\n');

            fname = fullfile(obj.topasFolder,obj.infilenames.geometry);
            matRad_cfg.dispInfo('Reading Geometry from %s\n',fname);
            world = fileread(fname);
            fprintf(fID,'%s\n',world);

        end

        function writeScorers(obj,fID)
            %TODO: Insert documentation
            matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class

            obj.MCparam.outputType = obj.scorer.outputType;

            % write dose to medium scorer
            if obj.scorer.doseToMedium
                fname = fullfile(obj.topasFolder,filesep,obj.scorerFolder,filesep,obj.infilenames.Scorer_doseToMedium);
                matRad_cfg.dispDebug('Reading doseToMedium scorer from %s\n',fname);
                scorerName = fileread(fname);
                fprintf(fID,'\n%s\n\n',scorerName);

                % Update MCparam.tallies with processed scorer
                obj.MCparam.tallies = [obj.MCparam.tallies,{'physicalDose'}];
            end

            % write RBE scorer
            if obj.scorer.RBE
                for i = 1:length(obj.scorer.RBE_model)
                    switch obj.radiationMode
                        case 'protons'
                            % Process available varRBE models for protons

                            if ~isempty(strfind(lower(obj.scorer.RBE_model{i}),'mcn'))
                                fname = fullfile(obj.topasFolder,filesep,obj.scorerFolder,filesep,obj.infilenames.Scorer_RBE_MCN);
                            elseif ~isempty(strfind(lower(obj.scorer.RBE_model{i}),'wed'))
                                fname = fullfile(obj.topasFolder,filesep,obj.scorerFolder,filesep,obj.infilenames.Scorer_RBE_WED);
                            else
                                matRad_cfg.dispError(['Model ',obj.scorer.RBE_model{i},' not implemented for ',obj.radiationMode]);
                            end
                        case {'carbon','helium'}
                            % Process available varRBE models for carbon and helium
                            if ~isempty(strfind(lower(obj.scorer.RBE_model{i}),'libamtrack'))
                                fname = fullfile(obj.topasFolder,filesep,obj.scorerFolder,filesep,obj.infilenames.Scorer_RBE_libamtrack);
                            elseif ~isempty(strfind(lower(obj.scorer.RBE_model{i}),'lem'))
                                fname = fullfile(obj.topasFolder,filesep,obj.scorerFolder,filesep,obj.infilenames.Scorer_RBE_LEM1);
                            else
                                matRad_cfg.dispError(['Model ',obj.scorer.RBE_model{i},' not implemented for ',obj.radiationMode]);
                            end
                        otherwise
                            % Throw error in case an invalid radiationMode has been selected
                            matRad_cfg.dispError(['Model ',obj.scorer.RBE_model{i},' not implemented for ',obj.radiationMode]);
                    end

                    % Read appropriate scorer from file and write to config file
                    matRad_cfg.dispDebug('Reading RBE Scorer from %s\n',fname);
                    scorerName = fileread(fname);
                    fprintf(fID,'\n%s\n\n',scorerName);
                end

                % Begin writing biological scorer components: cell lines
                switch obj.radiationMode
                    case 'protons'
                        fprintf(fID,'\n### Biological Parameters ###\n');
                        fprintf(fID,'sv:Sc/CellLines = 1 "CellLineGeneric"\n');
                        fprintf(fID,'d:Sc/CellLineGeneric/Alphax 		= Sc/AlphaX /Gy\n');
                        fprintf(fID,'d:Sc/CellLineGeneric/Betax 		= Sc/BetaX /Gy2\n');
                        fprintf(fID,'d:Sc/CellLineGeneric/AlphaBetaRatiox 	= Sc/AlphaBetaX Gy\n\n');
                    case {'carbon','helium'}
                        fprintf(fID,'\n### Biological Parameters ###\n');
                        fprintf(fID,'sv:Sc/CellLines = 1 "CellGeneric_abR2"\n');
                        fprintf(fID,'d:Sc/CellGeneric_abR2/Alphax = Sc/AlphaX /Gy\n');
                        fprintf(fID,'d:Sc/CellGeneric_abR2/Betax = Sc/BetaX /Gy2\n\n');
                        % fprintf(fID,'d:Sc/CellGeneric_abR2/AlphaBetaRatiox 	= Sc/AlphaBetaX Gy\n');
                    otherwise
                        matRad_cfg.dispError([obj.radiationMode ' not implemented']);
                end

                % write biological scorer components: dose parameters
                matRad_cfg.dispDebug('Writing Biologial Scorer components.\n');
                fprintf(fID,'d:Sc/PrescribedDose = %.4f Gy\n',obj.bioParameters.PrescribedDose);
                fprintf(fID,'b:Sc/SimultaneousExposure = %s\n',obj.bioParameters.SimultaneousExposure);
                fprintf(fID,'d:Sc/AlphaX = %.4f /Gy\n',obj.bioParameters.AlphaX);
                fprintf(fID,'d:Sc/BetaX = %.4f /Gy2\n',obj.bioParameters.BetaX);
                fprintf(fID,'d:Sc/AlphaBetaX = %.4f Gy\n',obj.bioParameters.AlphaX/obj.bioParameters.BetaX);

                % Update MCparam.tallies with processed scorer
                for i = 1:length(obj.scorer.RBE_model)
                    obj.MCparam.tallies = [obj.MCparam.tallies,{['alpha_' obj.scorer.RBE_model{i}],['beta_' obj.scorer.RBE_model{i}]}];
                end
            end

            % Write share sub-scorer
            if obj.scorer.sharedSubscorers && obj.scorer.RBE
                % Select appropriate scorer from selected flags
                scorerNames = {'Alpha','Beta'};
                if any(cellfun(@(teststr) ~isempty(strfind(lower(teststr),'mcn')), obj.scorer.RBE_model))
                    obj.scorer.LET = true;
                    obj.scorer.doseToWater = true;
                    scorerPrefix = 'McNamara';
                elseif any(cellfun(@(teststr) ~isempty(strfind(lower(teststr),'wed')), obj.scorer.RBE_model))
                    obj.scorer.LET = true;
                    obj.scorer.doseToWater = true;
                    scorerPrefix = 'Wedenberg';
                elseif any(cellfun(@(teststr) ~isempty(strfind(lower(teststr),'lem')), obj.scorer.RBE_model)) || any(cellfun(@(teststr) ~isempty(strfind(lower(teststr),'libamtrack')), obj.scorer.RBE_model))
                    obj.scorer.doseToWater = true;
                    scorerPrefix = 'tabulated';
                end

                % Write subscorer to config files
                for s = 1:length(scorerNames)
                    if strcmp(obj.radiationMode,'protons')
                        fprintf(fID,'s:Sc/%s%s/ReferencedSubScorer_LET      = "ProtonLET"\n',scorerPrefix,scorerNames{s});
                    end
                    fprintf(fID,'s:Sc/%s%s/ReferencedSubScorer_Dose     = "Tally_DoseToWater"\n',scorerPrefix,scorerNames{s});
                end
            end

            % write dose to water scorer from file
            if obj.scorer.doseToWater
                fname = fullfile(obj.topasFolder,filesep,obj.scorerFolder,filesep,obj.infilenames.Scorer_doseToWater);
                matRad_cfg.dispDebug('Reading doseToWater scorer from %s\n',fname);
                scorerName = fileread(fname);
                fprintf(fID,'\n%s\n\n',scorerName);

                % Update MCparam.tallies with processed scorer
                obj.MCparam.tallies = [obj.MCparam.tallies,{'doseToWater'}];
            end

            % write LET scorer from file
            if obj.scorer.LET
                if strcmp(obj.radiationMode,'protons')
                    fname = fullfile(obj.topasFolder,filesep,obj.scorerFolder,filesep,obj.infilenames.Scorer_LET);
                    matRad_cfg.dispDebug('Reading LET Scorer from %s\n',fname);
                    scorerName = fileread(fname);
                    fprintf(fID,'\n%s\n\n',scorerName);

                    % Update MCparam.tallies with processed scorer
                    obj.MCparam.tallies = [obj.MCparam.tallies,{'LET'}];
                else
                    matRad_cfg.dispError('LET in TOPAS only for protons!\n');
                end
            end

            % write volume scorer from file
            if obj.scorer.volume
                fileList = dir(fullfile(obj.topasFolder,filesep,obj.scorerFolder,filesep,'TOPAS_scorer_volume_*.in'));
                for fileIx=1:length(fileList)
                    fname = fullfile(obj.topasFolder,fileList(fileIx).name);
                    matRad_cfg.dispDebug('Reading Volume Scorer from %s\n',fname);
                    scorerName = fileread(fname);
                    fprintf(fID,'\n%s\n\n',scorerName);

                    tallyLabel = regexprep(fileList(fileIx).name,'TOPAS_scorer_volume_','');
                    tallyLabel = regexprep(tallyLabel,'.txt.in','');

                    % Update MCparam.tallies with processed scorer
                    obj.MCparam.tallies = [obj.MCparam.tallies,{tallyLabel}];
                end
            end

            % write surface track count from file
            if obj.scorer.surfaceTrackCount
                fname = fullfile(obj.topasFolder,filesep,obj.scorerFolder,filesep,obj.infilenames.Scorer_surfaceTrackCount);
                matRad_cfg.dispDebug('Reading surface scorer from %s\n',fname);
                scorerName = fileread(fname);
                fprintf(fID,'\n%s\n\n',scorerName);

                % Update MCparam.tallies with processed scorer
                obj.MCparam.tallies = [obj.MCparam.tallies,{'IC'}];
            end


            % Write timefeature-splitting in case of dij calculation
            if obj.scorer.calcDij
                tallyName = cell(1,0);
                if obj.scorer.RBE
                    if any(cellfun(@(teststr) ~isempty(strfind(lower(teststr),'mcn')), obj.MCparam.RBE_models))
                        tallyName{end+1} = 'McNamaraAlpha';
                        tallyName{end+1} = 'McNamaraBeta';
                    end
                    if any(cellfun(@(teststr) ~isempty(strfind(lower(teststr),'wed')), obj.MCparam.RBE_models))
                        tallyName{end+1} = 'WedenbergAlpha';
                        tallyName{end+1} = 'WedenbergBeta';
                    end
                    if any(cellfun(@(teststr) ~isempty(strfind(lower(teststr),'libamtrack')), obj.MCparam.RBE_models))
                        tallyName{end+1} = 'tabulatedAlpha';
                        tallyName{end+1} = 'tabulatedBeta';
                    end
                    if any(cellfun(@(teststr) ~isempty(strfind(lower(teststr),'lem')), obj.MCparam.RBE_models))
                        tallyName{end+1} = 'tabulatedAlpha';
                        tallyName{end+1} = 'tabulatedBeta';
                    end
                end
                if obj.scorer.LET
                    tallyName{end+1} = 'ProtonLET';
                end
                if obj.scorer.surfaceTrackCount
                    tallyName{end+1} = 'IC';
                end
                if obj.scorer.doseToMedium
                    tallyName{end+1} = 'Patient/Tally_DoseToMedium';
                end
                if obj.scorer.doseToMedium
                    tallyName{end+1} = 'Tally_DoseToWater';
                end

                % We should discuss here if that's something that has to be available for photons as well, turned off for now
                if ~strcmp(obj.radiationMode,'photons')
                    fprintf(fID,'#-- Time feature splitting for dij calculation\n');

                    for i = 1:length(tallyName)
                        fprintf(fID,['s:Sc/' tallyName{i} '/SplitByTimeFeature = "ImageName"\n']);
                    end
                end
            end
        end

        function writeStfFields(obj,ct,stf,w,baseData)
            %TODO: Insert documentation
            matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class

            isPhoton = false;

            switch obj.radiationMode
                case 'photons'
                    % if photons
                    isPhoton = true;
                    if any(ismember(obj.beamProfile,{'biGaussian','simple'}))
                        matRad_cfg.dispWarning('beamProfile "%s" not available for photons, switching to "%s" as default.',obj.beamProfile,obj.defaultPhotonBeamProfile);
                        obj.beamProfile = obj.defaultPhotonBeamProfile;
                    end

                otherwise
                    % if particles
                    if ~any(ismember(obj.beamProfile,{'biGaussian','simple'}))
                        matRad_cfg.dispWarning('beamProfile "%s" not available for particles, switching to "%s" as default.',obj.beamProfile,obj.defaultParticleBeamProfile);
                        obj.beamProfile = obj.defaultParticleBeamProfile;
                    end
            end

            %Bookkeeping
            obj.MCparam.nbFields = length(stf);

            %Sanity check
            if numel(w) ~= sum([stf(:).totalNumOfBixels])
                matRad_cfg.dispError('Given number of weights (#%d) doesn''t match bixel count in stf (#%d)',numel(w), sum([stf(:).totalNumOfBixels]));
            end

            nParticlesTotalBixel = round(obj.numParticlesPerHistory * w);
            nParticlesTotal = sum(nParticlesTotalBixel);
            maxParticlesBixel = obj.numParticlesPerHistory * max(w(:));
            minParticlesBixel = round(max([obj.minRelWeight*maxParticlesBixel,1]));

            switch obj.modeHistories
                case 'num'
                    obj.fracHistories = obj.numHistoriesDirect ./ sum(nParticlesTotalBixel);
                case 'frac'
                    obj.numHistoriesDirect = sum(nParticlesTotalBixel);
                otherwise
                    matRad_cfg.dispError('Invalid history setting!');
            end

            nParticlesTotal = 0;

            %Preread beam setup
            switch obj.beamProfile
                case 'biGaussian'
                    fname = fullfile(obj.topasFolder,obj.infilenames.beam_biGaussian);
                    TOPAS_beamSetup = fileread(fname);
                    matRad_cfg.dispInfo('Reading ''%s'' Beam Characteristics from ''%s''\n',obj.beamProfile,fname);

                case 'simple'
                    fname = fullfile(obj.topasFolder,obj.infilenames.beam_generic);
                    TOPAS_beamSetup = fileread(fname);
                    matRad_cfg.dispInfo('Reading ''%s'' Beam Characteristics from ''%s''\n',obj.beamProfile,fname);

                case 'phasespace'
                    fname = fullfile(obj.topasFolder,obj.infilenames.beam_phasespace);
                    TOPAS_beamSetup = fileread(fname);
                    obj.pencilBeamScanning = 0 ;
                    matRad_cfg.dispInfo('Reading ''%s'' Beam Characteristics from ''%s''\n',obj.beamProfile,fname);

                case 'virtualGaussian'
                    fname = fullfile(obj.topasFolder,obj.infilenames.beam_virtualGaussian);
                    TOPAS_beamSetup = fileread(fname);
                    matRad_cfg.dispInfo('Reading ''%s'' Beam Characteristics from ''%s''\n',obj.beamProfile,fname);

                case 'uniform'
                    fname = fullfile(obj.topasFolder,obj.infilenames.beam_uniform);
                    TOPAS_beamSetup = fileread(fname);
                    matRad_cfg.dispInfo('Reading ''%s'' Beam Characteristics from ''%s''\n',obj.beamProfile,fname);

                otherwise
                    matRad_cfg.dispError('Beam Type ''%s'' not supported for photons',obj.beamProfile);

            end

            % Set variables for loop over beams
            nBeamParticlesTotal = zeros(1,length(stf));
            currentBixel = 1;
            bixelNotMeetingParticleQuota = 0;
            historyCount = zeros(1,length(stf));

            for beamIx = 1:length(stf)

                SAD = stf(beamIx).SAD;

                if isPhoton
                    nozzleToAxisDistance = SAD;
                    sourceToNozzleDistance = 0;
                else
                    nozzleToAxisDistance = baseData.nozzleToIso;
                    sourceToNozzleDistance = SAD - nozzleToAxisDistance;

                    %Selection of base data given the energies and focusIndex
                    if obj.useOrigBaseData
                        [~,ixTmp,~] = intersect([ baseData.machine.data.energy], [stf.ray.energy]);
                        for i = 1:length(ixTmp)
                            selectedData(i) =  baseData.machine.data(ixTmp(i));
                        end
                        energies = [selectedData.energy];
                    else
                        selectedData = [];
                        focusIndex = baseData.selectedFocus(baseData.energyIndex);

                        scalarFields = ["NominalEnergy","EnergySpread","MeanEnergy"];                                                   
                        
                        for i = 1:numel(focusIndex)
                            for field = scalarFields
                                baseData.monteCarloData(i).(field) = ones(1,max(focusIndex))*baseData.monteCarloData(i).(field);
                            end
                            selectedData = [selectedData, structfun(@(x) x(focusIndex(i)),baseData.monteCarloData(i),'UniformOutput',false)];
                        end
                        energies = [selectedData.NominalEnergy];
                    end

                    %Get Range Shifters in field if present
                    allRays = [stf(beamIx).ray];
                    raShis = [allRays.rangeShifter];
                    [~,ix] =  unique(cell2mat(squeeze(struct2cell(raShis))'),'rows');

                    raShis = raShis(ix);
                    ix = [raShis.ID] == 0;
                    raShis = raShis(~ix);

                    %Convert ID into readable string
                    for r = 1:numel(raShis)
                        if isnumeric(raShis(r).ID)
                            raShis(r).topasID = ['RangeShifter' num2str(raShis(r).ID)];
                        else
                            raShis(r).topasID = ['RangeShifter' raShis(r).ID];
                        end
                    end
                end

                %get beamlet properties for each bixel in the stf and write it into dataTOPAS
                cutNumOfBixel = 0;

                % Clear dataTOPAS from the previous beam
                dataTOPAS = [];

                %Loop over rays and then over spots on ray
                for rayIx = 1:stf(beamIx).numOfRays
                    for bixelIx = 1:stf(beamIx).numOfBixelsPerRay(rayIx)

                        nCurrentParticles = nParticlesTotalBixel(currentBixel);

                        % check whether there are (enough) particles for beam delivery
                        if (nCurrentParticles>minParticlesBixel)

                            %                             collectBixelIdx(end+1) = bixelIx;
                            cutNumOfBixel = cutNumOfBixel + 1;
                            bixelEnergy = stf(beamIx).ray(rayIx).energy(bixelIx);

                            dataTOPAS(cutNumOfBixel).posX = stf(beamIx).ray(rayIx).rayPos_bev(3);
                            dataTOPAS(cutNumOfBixel).posY = stf(beamIx).ray(rayIx).rayPos_bev(1);

                            dataTOPAS(cutNumOfBixel).current = uint32(obj.fracHistories * nCurrentParticles / obj.numOfRuns);

                            if obj.pencilBeamScanning
                                % angleX corresponds to the rotation around the X axis necessary to move the spot in the Y direction
                                % angleY corresponds to the rotation around the Y' axis necessary to move the spot in the X direction
                                % note that Y' corresponds to the Y axis after the rotation of angleX around X axis
                                % note that Y translates to -Y for TOPAS
                                dataTOPAS(cutNumOfBixel).angleX = atan(dataTOPAS(cutNumOfBixel).posY / SAD);
                                dataTOPAS(cutNumOfBixel).angleY = atan(-dataTOPAS(cutNumOfBixel).posX ./ (SAD ./ cos(dataTOPAS(cutNumOfBixel).angleX)));
                                % Translate posX and posY to patient coordinates
                                dataTOPAS(cutNumOfBixel).posX = (dataTOPAS(cutNumOfBixel).posX / SAD)*(SAD-nozzleToAxisDistance);
                                dataTOPAS(cutNumOfBixel).posY = (dataTOPAS(cutNumOfBixel).posY / SAD)*(SAD-nozzleToAxisDistance);
                            end

                            switch obj.radiationMode
                                case {'protons','carbon','helium'}
                                    [~,ixTmp,~] = intersect(energies, bixelEnergy);
                                    if obj.useOrigBaseData
                                        dataTOPAS(cutNumOfBixel).energy = selectedData(ixTmp).energy;
                                        dataTOPAS(cutNumOfBixel).focusFWHM = selectedData(ixTmp).initFocus.SisFWHMAtIso(stf(beamIx).ray(rayIx).focusIx(bixelIx));

                                    else
                                        dataTOPAS(cutNumOfBixel).energy = selectedData(ixTmp).MeanEnergy;
                                        dataTOPAS(cutNumOfBixel).nominalEnergy = selectedData(ixTmp).NominalEnergy;
                                        dataTOPAS(cutNumOfBixel).energySpread = selectedData(ixTmp).EnergySpread;
                                        dataTOPAS(cutNumOfBixel).spotSizeX = selectedData(ixTmp).SpotSize1x;
                                        dataTOPAS(cutNumOfBixel).divergenceX = selectedData(ixTmp).Divergence1x;
                                        dataTOPAS(cutNumOfBixel).correlationX = selectedData(ixTmp).Correlation1x;
                                        dataTOPAS(cutNumOfBixel).spotSizeY = selectedData(ixTmp).SpotSize1y;
                                        dataTOPAS(cutNumOfBixel).divergenceY = selectedData(ixTmp).Divergence1y;
                                        dataTOPAS(cutNumOfBixel).correlationY = selectedData(ixTmp).Correlation1y;
                                        dataTOPAS(cutNumOfBixel).focusFWHM = selectedData(ixTmp).FWHMatIso;
                                    end
                                case 'photons'
                                    dataTOPAS(cutNumOfBixel).energy = bixelEnergy;
                                    dataTOPAS(cutNumOfBixel).energySpread = 0;
                            end

                            if obj.scorer.calcDij
                                % remember beam and bixel number
                                dataTOPAS(cutNumOfBixel).beam           = beamIx;
                                dataTOPAS(cutNumOfBixel).ray            = rayIx;
                                dataTOPAS(cutNumOfBixel).bixel          = bixelIx;
                                dataTOPAS(cutNumOfBixel).totalBixel     = currentBixel;
                            end

                            %Add RangeShifterState
                            if exist('raShis','var') && ~isempty(raShis)
                                raShiOut = zeros(1,length(raShis));
                                for r = 1:length(raShis)
                                    if stf(beamIx).ray(rayIx).rangeShifter(bixelIx).ID == raShis(r).ID
                                        raShiOut(r) = 0; %Range shifter is in beam path
                                    else
                                        raShiOut(r) = 1; %Range shifter is out of beam path / not used
                                    end
                                end
                                dataTOPAS(cutNumOfBixel).raShiOut = raShiOut;
                            end

                            nBeamParticlesTotal(beamIx) = nBeamParticlesTotal(beamIx) + nCurrentParticles;


                        end

                        currentBixel = currentBixel + 1;

                    end
                end

                bixelNotMeetingParticleQuota = bixelNotMeetingParticleQuota + (stf(beamIx).totalNumOfBixels-cutNumOfBixel);

                % discard data if the current has unphysical values
                idx = find([dataTOPAS.current] < 1);
                dataTOPAS(idx) = [];

                % Safety check for empty beam (not allowed)
                if isempty(dataTOPAS)
                    matRad_cfg.dispError('dataTOPAS of beam %i is empty.',beamIx);
                else
                    cutNumOfBixel = length(dataTOPAS(:));
                end

                % Sort dataTOPAS according to energy
                if length(dataTOPAS)>1 && ~issorted([dataTOPAS(:).energy])
                    [~,ixSorted] = sort([dataTOPAS(:).energy]);
                    dataTOPAS = dataTOPAS(ixSorted);
                end

                % Save adjusted beam histories
                historyCount(beamIx) = uint32(obj.fracHistories * nBeamParticlesTotal(beamIx) / obj.numOfRuns);

                if historyCount(beamIx) < cutNumOfBixel || cutNumOfBixel == 0
                    matRad_cfg.dispError('Insufficient number of histories!')
                end

                % Check if current has the set amount of histories
                % If needed, adjust current to actual histories (by adding/subtracting from random rays)
                while sum([dataTOPAS(:).current]) ~= historyCount(beamIx)
                    diff = sum([dataTOPAS.current]) - sum(historyCount(beamIx));
                    if matRad_cfg.isMatlab
                        [~,~,R] = histcounts(rand(abs(diff),1),cumsum([0;double(transpose([dataTOPAS(:).current]))./double(sum([dataTOPAS(:).current]))]));
                    else
                        [~,R] = histc(rand(abs(diff),1),cumsum([0;double(transpose([dataTOPAS(:).current]))./double(sum([dataTOPAS(:).current]))]));
                    end
                    idx = 1:length(dataTOPAS);
                    randIx = idx(R);

                    newCurr = num2cell(arrayfun(@plus,double([dataTOPAS(randIx).current]),-1*sign(diff)*ones(1,abs(diff))),1);
                    [dataTOPAS(randIx).current] = newCurr{:};
                end

                % Previous histories were set per run
                historyCount(beamIx) = historyCount(beamIx) * obj.numOfRuns;

                % Write TOPAS data base file
                if isfield(ct,'currCtScen')
                    % 4D case
                    fieldSetupFileName = sprintf('beamSetup_%s_field%d_ct%d.txt',obj.label,beamIx,ct.currCtScen);
                    fileID = fopen(fullfile(obj.workingDir,fieldSetupFileName),'w');
                    obj.writeFieldHeader(fileID,ct.currCtScen);
                else
                    fieldSetupFileName = sprintf('beamSetup_%s_field%d.txt',obj.label,beamIx);
                    fileID = fopen(fullfile(obj.workingDir,fieldSetupFileName),'w');
                    obj.writeFieldHeader(fileID);
                end

                % NozzleAxialDistance
                if isPhoton 
                    fprintf(fileID,'d:Ge/Nozzle/TransZ = -%f mm\n', stf(beamIx).SCD+40); %Phasespace hardcorded infront of MLC at SSD 46 cm
                else
                    fprintf(fileID,'d:Ge/Nozzle/TransZ = -%f mm\n', nozzleToAxisDistance);
                end

                if obj.pencilBeamScanning
                    fprintf(fileID,'d:Ge/Nozzle/RotX = Tf/Beam/AngleX/Value rad\n');
                    fprintf(fileID,'d:Ge/Nozzle/RotY = Tf/Beam/AngleY/Value rad\n');
                    fprintf(fileID,'d:Ge/Nozzle/RotZ = 0.0 rad\n\n');
                end

                %Write modality specific info
                switch stf(beamIx).radiationMode
                    case 'protons'
                        fprintf(fileID,'s:Sim/ParticleName = "proton"\n');
                        fprintf(fileID,'u:Sim/ParticleMass = 1.0\n');

                        particleA = 1;
                        % particleZ = 1;

                        modules = obj.modules_protons;

                    case 'carbon'
                        fprintf(fileID,'s:Sim/ParticleName = "GenericIon(6,12)"\n');
                        fprintf(fileID,'u:Sim/ParticleMass = 12.0\n');

                        particleA = 12;
                        % particleZ = 6;

                        modules = obj.modules_GenericIon;

                    case 'helium'
                        fprintf(fileID,'s:Sim/ParticleName = "GenericIon(2,4)"\n');
                        fprintf(fileID,'u:Sim/ParticleMass = 4.0\n');

                        particleA = 4;
                        % particleZ = 2;

                        modules = obj.modules_GenericIon;

                    case 'photons'
                        fprintf(fileID,'s:Sim/ParticleName = "gamma"\n');
                        fprintf(fileID,'u:Sim/ParticleMass = 0\n');

                        particleA = 0;
                        % particleZ = 0;

                        modules = obj.modules_photons;

                    otherwise
                        matRad_cfg.dispError('Invalid radiation mode %s!',stf.radiationMode)
                end

                if obj.pencilBeamScanning
                    % Write couch and gantry angles
                    fprintf(fileID,'d:Sim/GantryAngle = %f deg\n',stf(beamIx).gantryAngle);
                    fprintf(fileID,'d:Sim/CouchAngle = %f deg\n',stf(beamIx).couchAngle);

                    % Write time feature (TOPAS uses time features to loop through bixels)
                    fprintf(fileID,'d:Tf/TimelineStart = 0. ms\n');
                    fprintf(fileID,'d:Tf/TimelineEnd = %i ms\n', 10 * cutNumOfBixel);
                    fprintf(fileID,'i:Tf/NumberOfSequentialTimes = %i\n', cutNumOfBixel);
                    fprintf(fileID,'dv:Tf/Beam/Spot/Times = %i ', cutNumOfBixel);
                    fprintf(fileID,'%i ',linspace(10,cutNumOfBixel*10,cutNumOfBixel));
                    fprintf(fileID,' ms\n');
                    %fprintf(fileID,'uv:Tf/Beam/Spot/Values = %i %s\n',cutNumOfBixel,num2str(collectBixelIdx));

                    % Write energySpectrum if available and flag is set
                    if ~isPhoton && isfield(baseData.machine.data,'energySpectrum') && obj.useEnergySpectrum
                        matRad_cfg.dispInfo('Beam energy spectrum available\n');
                        energySpectrum = [baseData.machine.data(:).energySpectrum];
                        nbSpectrumPoints = length(energySpectrum(1).energy_MeVpN);

                        % Get energy indices of the current energies in the baseData
                        [~,energyIx] = ismember([dataTOPAS.nominalEnergy],[baseData.machine.data.energy]);

                        fprintf(fileID,'s:So/PencilBeam/BeamEnergySpectrumType = "Continuous"\n');
                        fprintf(fileID,'dv:So/PencilBeam/BeamEnergySpectrumValues = %d %s MeV\n',nbSpectrumPoints,strtrim(sprintf('Tf/Beam/EnergySpectrum/Energy/Point%03d/Value ',1:nbSpectrumPoints)));
                        fprintf(fileID,'uv:So/PencilBeam/BeamEnergySpectrumWeights = %d %s\n',nbSpectrumPoints,strtrim(sprintf('Tf/Beam/EnergySpectrum/Weight/Point%03d/Value ',1:nbSpectrumPoints)));
                        points_energy = reshape([energySpectrum(energyIx).energy_MeVpN],[],length(energyIx));
                        points_weight = reshape([energySpectrum(energyIx).weight],[],length(energyIx));
                        for spectrumPoint=1:nbSpectrumPoints
                            fprintf(fileID,'s:Tf/Beam/EnergySpectrum/Energy/Point%03d/Function = "Step"\n',spectrumPoint);
                            fprintf(fileID,'dv:Tf/Beam/EnergySpectrum/Energy/Point%03d/Times = Tf/Beam/Spot/Times ms\n',spectrumPoint);
                            fprintf(fileID,'dv:Tf/Beam/EnergySpectrum/Energy/Point%03d/Values = %d %s MeV\n',spectrumPoint,cutNumOfBixel,strtrim(sprintf('%f ',particleA*points_energy(spectrumPoint,:))));
                            fprintf(fileID,'s:Tf/Beam/EnergySpectrum/Weight/Point%03d/Function = "Step"\n',spectrumPoint);
                            fprintf(fileID,'dv:Tf/Beam/EnergySpectrum/Weight/Point%03d/Times = Tf/Beam/Spot/Times ms\n',spectrumPoint);
                            fprintf(fileID,'uv:Tf/Beam/EnergySpectrum/Weight/Point%03d/Values = %d %s\n',spectrumPoint,cutNumOfBixel,strtrim(sprintf('%f ',points_weight(spectrumPoint,:))));
                        end
                    end

                    % Write amount of energies in plan
                    fprintf(fileID,'s:Tf/Beam/Energy/Function = "Step"\n');
                    fprintf(fileID,'dv:Tf/Beam/Energy/Times = Tf/Beam/Spot/Times ms\n');
                    fprintf(fileID,'dv:Tf/Beam/Energy/Values = %i ', cutNumOfBixel);

                    % Write actual energies
                    % WARNING: Transform total energy with atomic number
                    fprintf(fileID,'%f ',particleA*[dataTOPAS.energy]);
                    fprintf(fileID,' MeV\n');
                end

                % Write beam profile
                switch obj.beamProfile
                    case 'biGaussian'
                        % Write energy spread
                        fprintf(fileID,'s:Tf/Beam/EnergySpread/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/EnergySpread/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/EnergySpread/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,'%f ',[dataTOPAS.energySpread]);
                        fprintf(fileID,'\n');

                        % Write parameters for first dimension
                        fprintf(fileID,'s:Tf/Beam/SigmaX/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaX/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaX/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,'%f ',[dataTOPAS.spotSizeX]);
                        fprintf(fileID,' mm\n');
                        fprintf(fileID,'s:Tf/Beam/SigmaXPrime/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaXPrime/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/SigmaXPrime/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,'%f ',[dataTOPAS.divergenceX]);
                        fprintf(fileID,'\n');
                        fprintf(fileID,'s:Tf/Beam/CorrelationX/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/CorrelationX/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/CorrelationX/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,'%f ',[dataTOPAS.correlationX]);
                        fprintf(fileID,'\n');

                        % Write parameters for second dimension (profile is uniform)
                        fprintf(fileID,'s:Tf/Beam/SigmaY/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaY/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaY/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,'%f ',[dataTOPAS.spotSizeY]);
                        fprintf(fileID,' mm\n');
                        fprintf(fileID,'s:Tf/Beam/SigmaYPrime/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaYPrime/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/SigmaYPrime/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,'%f ',[dataTOPAS.divergenceY]);
                        fprintf(fileID,'\n');
                        fprintf(fileID,'s:Tf/Beam/CorrelationY/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/CorrelationY/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/CorrelationY/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,'%f ',[dataTOPAS.correlationY]);
                        fprintf(fileID,'\n');

                    case 'simple'
                        fprintf(fileID,'s:Tf/Beam/FocusFWHM/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/FocusFWHM/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'dv:Tf/Beam/FocusFWHM/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,'%f ',[dataTOPAS.focusFWHM]);
                        fprintf(fileID,' mm\n');

                    case 'virtualGaussian'
                        fprintf(fileID,'s:Tf/Beam/EnergySpread/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/EnergySpread/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/EnergySpread/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.energySpread]));
                        fprintf(fileID,'\n');

                        if isfield([stf.ray], 'collimation')
                            % Use field width for now
                            fprintf(fileID,'d:So/PencilBeam/BeamPositionSpreadX = %d mm\n', stf(1).ray.collimation.fieldWidth);
                            fprintf(fileID,'d:So/PencilBeam/BeamPositionSpreadY = %d mm\n', stf(1).ray.collimation.fieldWidth);
                        else
                            % Set some default value
                            fprintf(fileID,'d:So/PencilBeam/BeamPositionSpreadX = %d mm\n', 30);
                            fprintf(fileID,'d:So/PencilBeam/BeamPositionSpreadY = %d mm\n', 30);
                        end

                    case 'uniform'
                        fprintf(fileID,'s:Tf/Beam/EnergySpread/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/EnergySpread/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/EnergySpread/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.energySpread]));
                        fprintf(fileID,'\n');

                        if  isfield([stf.ray],'collimation')
                            % Use field width for now
                            fprintf(fileID,'d:So/PencilBeam/BeamPositionCutoffX = %d mm\n', stf(1).ray.collimation.fieldWidth/2);
                            fprintf(fileID,'d:So/PencilBeam/BeamPositionCutoffY = %d mm\n', stf(1).ray.collimation.fieldWidth/2);
                        else
                            % Set some default value
                            fprintf(fileID,'d:So/PencilBeam/BeamPositionCutoffX = %d mm\n', 15);
                            fprintf(fileID,'d:So/PencilBeam/BeamPositionCutoffY = %d mm\n', 15);
                        end

                    case 'phasespace'

                        fprintf(fileID,'d:Sim/GantryAngle = %f deg\n',stf(beamIx).gantryAngle); %just one beam angle for now
                        fprintf(fileID,'d:Sim/CouchAngle = %f deg\n',stf(beamIx).couchAngle);
                        % Here the phasespace file is loaded and referenced in the beamSetup file                      
                        if strcmp(obj.externalCalculation,'write')
                            matRad_cfg.dispWarning(['External calculation and phaseSpace selected, manually place ' obj.infilenames.phaseSpaceSourcePhotons '.header and ' obj.infilenames.phaseSpaceSourcePhotons  '.phsp into your simulation directory.']);
                        else
                            if length(dir([obj.topasFolder filesep 'beamSetup' filesep 'phasespace' filesep obj.infilenames.phaseSpaceSourcePhotons '*'])) < 2
                                matRad_cfg.dispError([phaseSpaceFileName ' header or phsp file could not be found in beamSetup/phasespace folder.']);
                            end
                        end
                        %phasespaceStr = ['..' filesep 'beamSetup' filesep 'phasespace' filesep phaseSpaceFileName];
                        %&phasespaceStr =  replace(phasespaceStr, '\', '/');
                        fprintf(fileID,'s:So/Phasespace/PhaseSpaceFileName = "%s"\n', obj.infilenames.phaseSpaceSourcePhotons );

                end

                % Write spot angles
                if obj.pencilBeamScanning
                    fprintf(fileID,'s:Tf/Beam/AngleX/Function = "Step"\n');
                    fprintf(fileID,'dv:Tf/Beam/AngleX/Times = Tf/Beam/Spot/Times ms\n');
                    fprintf(fileID,'dv:Tf/Beam/AngleX/Values = %i ', cutNumOfBixel);
                    fprintf(fileID,'%f ',[dataTOPAS.angleX]);
                    fprintf(fileID,' rad\n');
                    fprintf(fileID,'s:Tf/Beam/AngleY/Function = "Step"\n');
                    fprintf(fileID,'dv:Tf/Beam/AngleY/Times = Tf/Beam/Spot/Times ms\n');
                    fprintf(fileID,'dv:Tf/Beam/AngleY/Values = %i ', cutNumOfBixel);
                    fprintf(fileID,'%f ',[dataTOPAS.angleY]);
                    fprintf(fileID,' rad\n');


                    % Write spot positions
                    fprintf(fileID,'s:Tf/Beam/PosX/Function = "Step"\n');
                    fprintf(fileID,'dv:Tf/Beam/PosX/Times = Tf/Beam/Spot/Times ms\n');
                    fprintf(fileID,'dv:Tf/Beam/PosX/Values = %i ', cutNumOfBixel);
                    fprintf(fileID,'%f ',[dataTOPAS.posX]);
                    fprintf(fileID,' mm\n');
                    fprintf(fileID,'s:Tf/Beam/PosY/Function = "Step"\n');
                    fprintf(fileID,'dv:Tf/Beam/PosY/Times = Tf/Beam/Spot/Times ms\n');
                    fprintf(fileID,'dv:Tf/Beam/PosY/Values = %i ', cutNumOfBixel);
                    fprintf(fileID,'%f ',[dataTOPAS.posY]);
                    fprintf(fileID,' mm\n');

                    % Write spot current (translates to the amount of particles in a spot)
                    fprintf(fileID,'s:Tf/Beam/Current/Function = "Step"\n');
                    fprintf(fileID,'dv:Tf/Beam/Current/Times = Tf/Beam/Spot/Times ms\n');
                    fprintf(fileID,'iv:Tf/Beam/Current/Values = %i ', cutNumOfBixel);
                    fprintf(fileID,'%i ',[dataTOPAS.current]);
                    fprintf(fileID,'\n\n');

                    % Range shifter in/out
                    if ~isPhoton && ~isempty(raShis)
                        fprintf(fileID,'#Range Shifter States:\n');
                        for r = 1:numel(raShis)
                            fprintf(fileID,'s:Tf/Beam/%sOut/Function = "Step"\n',raShis(r).topasID);
                            fprintf(fileID,'dv:Tf/Beam/%sOut/Times = Tf/Beam/Spot/Times ms\n',raShis(r).topasID);
                            fprintf(fileID,'uv:Tf/Beam/%sOut/Values = %i ', raShis(r).topasID, cutNumOfBixel);
                            fprintf(fileID,'%f ',[dataTOPAS.raShiOut]);
                            fprintf(fileID,'\n\n');
                        end

                        % Range Shifter Definition
                        for r = 1:numel(raShis)
                            obj.writeRangeShifter(fileID,raShis(r),sourceToNozzleDistance);
                        end
                    end


                end

                % Write previously beam profile
                fprintf(fileID,'%s\n',TOPAS_beamSetup);

                % Write MLC if available
                if isfield(stf(beamIx).ray, 'shapes')
                    SCD = stf(beamIx).SCD;
                    fname = fullfile(obj.topasFolder,obj.infilenames.beam_mlc);
                    TOPAS_mlcSetup = fileread(fname);
                    fprintf(fileID,'%s\n',TOPAS_mlcSetup);
                    [numOfLeaves,leafTimes]=size([stf(beamIx).ray.shapes(:).leftLeafPos]); %there are #numOfLeaves leaves and #leafTimes times/shapes
                    leftLeafPos = [stf(beamIx).ray.shapes(:).leftLeafPos]*SCD./SAD;
                    rightLeafPos = [stf(beamIx).ray.shapes(:).rightLeafPos]*SCD./SAD;
                    % Set MLC paramters as in TOPAS example file https://topas.readthedocs.io/en/latest/parameters/geometry/specialized.html#multi-leaf-collimator
                    fprintf(fileID,'d:Sim/Ge/MultiLeafCollimatorA/TransZ   = %f cm\n', 4);
                    fprintf(fileID,'d:Ge/MultiLeafCollimatorA/MaximumLeafOpen   = %f cm\n',15);
                    fprintf(fileID,'d:Ge/MultiLeafCollimatorA/Thickness         = %f cm\n',8);
                    fprintf(fileID,'d:Ge/MultiLeafCollimatorA/Length            = %f  cm\n',15);
                    fprintf(fileID,'dv:Ge/MultiLeafCollimatorA/Widths           = %i ', numOfLeaves+2);
                    fprintf(fileID, '%f ', [200,  stf(1).ray.collimation.leafWidth*ones(1,numOfLeaves)*SCD./SAD , 200]);
                    fprintf(fileID,' mm \n');
                    fprintf(fileID,'dv:Ge/MultiLeafCollimatorA/XPlusLeavesOpen  = %i ',numOfLeaves+2);
                    for i = 0:numOfLeaves+1
                        fprintf( fileID,'Tf/LeafXPlus%i/Value ',i);
                    end
                    fprintf(fileID,'mm \n');
                    fprintf(fileID,'dv:Ge/MultiLeafCollimatorA/XMinusLeavesOpen  = %i ',numOfLeaves+2);
                    for i = 0:numOfLeaves+1
                        fprintf( fileID,'Tf/LeafXMinus%i/Value ',i);
                    end
                    fprintf(fileID,'mm \n');

                    %initilization of time features
                    fprintf(fileID,'d:Tf/TimelineStart = 0 ms\n');
                    fprintf(fileID,'d:Tf/TimelineEnd = %f ms\n',leafTimes*10);
                    fprintf(fileID,'i:Tf/NumberOfSequentialTimes = %i \n',leafTimes);

                    for i = 1:numOfLeaves
                        fprintf(fileID,'s:Tf/LeafXMinus%i/Function  = "Step"\n',i);
                        fprintf(fileID,'dv:Tf/LeafXMinus%i/Times =  %i ', i,leafTimes);
                        fprintf(fileID,'%i ', [1:leafTimes]*10);
                        fprintf(fileID,' ms\n');
                        fprintf(fileID,'dv:Tf/LeafXMinus%i/Values = %i ', i,leafTimes);
                        fprintf(fileID,'%f ', leftLeafPos(i,:));
                        fprintf(fileID,' mm\n\n');

                        fprintf(fileID,'s:Tf/LeafXPlus%i/Function  = "Step"\n',i);
                        fprintf(fileID,'dv:Tf/LeafXPlus%i/Times =  %i ',i,leafTimes);
                        fprintf(fileID,'%i ',[1:leafTimes]*10);
                        fprintf(fileID,' ms\n');
                        fprintf(fileID,'dv:Tf/LeafXPlus%i/Values = %i ', i,leafTimes);
                        fprintf(fileID,'%f ', rightLeafPos(i,:));
                        fprintf(fileID,' mm\n\n');
                    end
                    %Add aditional Leaf at the top and bottom to catch
                    %scattering
                    for i = [0,numOfLeaves+1]
                        fprintf(fileID,'s:Tf/LeafXMinus%i/Function  = "Step"\n',i);
                        fprintf(fileID,'dv:Tf/LeafXMinus%i/Times =  %i ', i,leafTimes);
                        fprintf(fileID,'%i ',[1:leafTimes]*10);
                        fprintf(fileID,' ms\n');
                        fprintf(fileID,'dv:Tf/LeafXMinus%i/Values = %i ', i,leafTimes);
                        fprintf(fileID,'%f ', zeros(size([1:leafTimes])));
                        fprintf(fileID,' mm\n\n');

                        fprintf(fileID,'s:Tf/LeafXPlus%i/Function  = "Step"\n',i);
                        fprintf(fileID,'dv:Tf/LeafXPlus%i/Times =  %i ',i,leafTimes);
                        fprintf(fileID,'%i ',[1:leafTimes]*10);
                        fprintf(fileID,' ms\n');
                        fprintf(fileID,'dv:Tf/LeafXPlus%i/Values = %i ', i,leafTimes);
                        fprintf(fileID,'%f ', zeros(size([1:leafTimes])));
                        fprintf(fileID,' mm\n\n');
                    end

                    fprintf(fileID, 's:Tf/Phasespace/NumberOfHistoriesInRun/Function  = "Step" \n');
                    fprintf(fileID, 'dv:Tf/Phasespace/NumberOfHistoriesInRun/Times = %i ', leafTimes);
                    fprintf(fileID,'%i ',[1:leafTimes]*10);
                    fprintf(fileID,' ms\n');
                    fprintf(fileID, 'iv:Tf/Phasespace/NumberOfHistoriesInRun/Values = %i ', leafTimes);
                    fprintf(fileID,'%i ',[dataTOPAS(:).current]);
                    fprintf(fileID,' \n');

                end

                % Translate patient according to beam isocenter
                fprintf(fileID,'d:Ge/Patient/TransX      = %f mm\n',0.5*ct.resolution.x*(ct.cubeDim(2)+1)-stf(beamIx).isoCenter(1));
                fprintf(fileID,'d:Ge/Patient/TransY      = %f mm\n',0.5*ct.resolution.y*(ct.cubeDim(1)+1)-stf(beamIx).isoCenter(2));
                fprintf(fileID,'d:Ge/Patient/TransZ      = %f mm\n',0.5*ct.resolution.z*(ct.cubeDim(3)+1)-stf(beamIx).isoCenter(3));
                fprintf(fileID,'d:Ge/Patient/RotX=0. deg\n');
                fprintf(fileID,'d:Ge/Patient/RotY=0. deg\n');
                fprintf(fileID,'d:Ge/Patient/RotZ=0. deg\n');

                % Load topas modules depending on the particle type
                fprintf(fileID,'\n# MODULES\n');
                moduleString = cellfun(@(s) sprintf('"%s"',s),modules,'UniformOutput',false);
                fprintf(fileID,'sv:Ph/Default/Modules = %d %s\n',length(modules),strjoin(moduleString,' '));

                fclose(fileID);
                % Write run scripts for TOPAS
                for runIx = 1:obj.numOfRuns
                    if isfield(ct,'currCtScen')
                        runFileName = sprintf('%s_field%d_ct%d_run%d.txt',obj.label,beamIx,ct.currCtScen,runIx);
                    else
                        runFileName = sprintf('%s_field%d_run%d.txt',obj.label,beamIx,runIx);
                    end
                    fileID = fopen(fullfile(obj.workingDir,runFileName),'w');

                    % Write header
                    if isfield(ct,'currCtScen')
                        obj.writeRunHeader(fileID,beamIx,runIx,ct.currCtScen);
                    else
                        obj.writeRunHeader(fileID,beamIx,runIx);
                    end

                    % Include path to beamSetup file
                    fprintf(fileID,'includeFile = ./%s\n',fieldSetupFileName);

                    % Write lines from scorer files
                    obj.writeScorers(fileID);

                    % Write dij-related config lines
                    % TODO: move this to github issue/todo -> We should discuss here if that's something that has to be available for photons as well
                    if ~strcmp(obj.radiationMode,'photons')
                        if obj.scorer.calcDij
                            fprintf(fileID,'\n');
                            fprintf(fileID,'#-- time feature splitting for dij calculation\n');
                            fprintf(fileID,'s:Tf/ImageName/Function = "Step"\n');
                            % create time feature scorer and save with original rays and bixel names
                            imageName = ['sv:Tf/ImageName/Values = ',num2str(cutNumOfBixel),cell2mat(strcat(strcat(' "ray',strsplit(num2str([dataTOPAS.ray]))),strcat('_bixel',strsplit(num2str([dataTOPAS.bixel])),'"')))];
                            fprintf(fileID,'%s\n',strjoin(imageName));
                            fprintf(fileID,'dv:Tf/ImageName/Times = Tf/Beam/Spot/Times ms\n');
                        end
                    end

                    fclose(fileID);
                end
            end

            if bixelNotMeetingParticleQuota ~= 0
                matRad_cfg.dispWarning([num2str(bixelNotMeetingParticleQuota) ' bixels were discarded due to particle threshold.'])
            end

            % Bookkeeping
            obj.MCparam.nbParticlesTotal = sum(nBeamParticlesTotal);
            obj.MCparam.nbHistoriesTotal = sum(historyCount);
            obj.MCparam.nbParticlesField = nBeamParticlesTotal;
            obj.MCparam.nbHistoriesField = historyCount;
        end

        function writePatient(obj,ct,pln)
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % matRad export CT RSP data for TOPAS simulation
            %
            % call
            %   obj.writePatient(ct, path, material)
            %
            % input
            %   ct:                 ct cube
            %   pln:                plan structure containing doseCalc classes
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % the image cube contains the indexing of materials
            % since at the moment TOPAS does not support ushort
            % the materials should have indexes between 0 and 32767
            % therefore, the maximum length of the vector is 32768

            matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class
            medium = obj.rsp_basematerial;
            if isequal(obj.arrayOrdering,'C')
                if matRad_cfg.logLevel > 2
                    matRad_cfg.dispInfo('Exporting cube in C ordering...\n')
                end
                permutation = [3 1 2];
            else
                if matRad_cfg.logLevel > 2
                    matRad_cfg.dispInfo('Exporting cube in FORTRAN ordering...\n')
                end
                permutation = [2 1 3];
            end

            % Bookkeeping
            obj.MCparam.imageCubeOrdering = obj.arrayOrdering;
            obj.MCparam.imageCubeConversionType = obj.materialConverter.mode;
            obj.MCparam.imageCubeFile = obj.outfilenames.patientCube;
            obj.MCparam.imageCubeDim = ct.cubeDim;
            obj.MCparam.imageVoxelDimension = ct.resolution;

            % Save filenames
            paramFile = obj.outfilenames.patientParam;
            dataFile = obj.outfilenames.patientCube;

            % Add ctScen number to filenames
            if isfield(ct,'currCtScen')
                ctScen = ct.currCtScen;
                paramFile = strsplit(paramFile,'.');
                paramFile = strjoin(paramFile,[num2str(ct.currCtScen) '.']);

                dataFile = strsplit(dataFile,'.');
                dataFile = strjoin(dataFile,[num2str(ct.currCtScen) '.']);
            else
                ctScen = 1;
            end

            % Open file to write in data
            outfile = fullfile(obj.workingDir, paramFile);
            matRad_cfg.dispInfo('Writing data to %s\n',outfile)
            fID = fopen(outfile,'w+');

            % Write material converter
            switch obj.materialConverter.mode
                case 'RSP' % Relative stopping power converter
                    rspHlut = matRad_loadHLUT(ct,obj.radiationMode);
                    min_HU = rspHlut(1,1);
                    max_HU = rspHlut(end,1);

                    huCube = int32(permute(ct.cubeHU{ctScen},permutation)); %  X,Y,Z ordering
                    huCube(huCube < min_HU) = min_HU;
                    huCube(huCube > max_HU) = max_HU;

                    unique_hu = unique(huCube(:));
                    unique_rsp = matRad_interp1(rspHlut(:,1),rspHlut(:,2),double(unique_hu));
                    fbase = fopen(['materialConverter/definedMaterials/' medium '.txt'],'r');
                    while ~feof(fbase)
                        strLine = fgets(fbase); %# read line by line
                        fprintf(fID,'%s',strLine);
                    end
                    fclose(fbase);

                    unique_materials = cell(1,length(unique_hu));
                    for ix=1:length(unique_hu)
                        unique_materials{ix} = strrep(['Material_HU_',num2str(unique_hu(ix))],'-','m');
                        fprintf(fID,'s:Ma/%s/BaseMaterial = "%s"\n',unique_materials{ix},medium);
                        fprintf(fID,'d:Ma/%s/Density = %f g/cm3\n',unique_materials{ix},unique_rsp(ix));
                    end

                    fprintf(fID,'s:Ge/Patient/Parent="World"\n');
                    fprintf(fID,'s:Ge/Patient/Type = "TsImageCube"\n');
                    fprintf(fID,'s:Ge/Patient/InputDirectory = "./"\n');
                    fprintf(fID,'s:Ge/Patient/InputFile = "%s"\n',dataFile);
                    fprintf(fID,'s:Ge/Patient/ImagingtoMaterialConverter = "MaterialTagNumber"\n');
                    fprintf(fID,'i:Ge/Patient/NumberOfVoxelsX = %d\n',ct.cubeDim(2));
                    fprintf(fID,'i:Ge/Patient/NumberOfVoxelsY = %d\n',ct.cubeDim(1));
                    fprintf(fID,'iv:Ge/Patient/NumberOfVoxelsZ = 1 %d\n',ct.cubeDim(3));
                    fprintf(fID,'d:Ge/Patient/VoxelSizeX       = %.3f mm\n',ct.resolution.x);
                    fprintf(fID,'d:Ge/Patient/VoxelSizeY       = %.3f mm\n',ct.resolution.y);
                    fprintf(fID,'dv:Ge/Patient/VoxelSizeZ       = 1 %.3f mm\n',ct.resolution.z);
                    fprintf(fID,'s:Ge/Patient/DataType  = "SHORT"\n');
                    fprintf(fID,'iv:Ge/Patient/MaterialTagNumbers = %d ',length(unique_hu));
                    fprintf(fID,num2str(unique_hu','%d '));
                    fprintf(fID,'\n');
                    fprintf(fID,'sv:Ge/Patient/MaterialNames = %d ',length(unique_hu));
                    fprintf(fID,'"%s"',strjoin(unique_materials,'" "'));
                    fprintf(fID,'\n');
                    fclose(fID);

                    % write data
                    fID = fopen(fullfile(obj.workingDir, dataFile),'w');
                    fwrite(fID,huCube,'short');
                    fclose(fID);
                    cube = huCube;


                case 'HUToWaterSchneider' % Schneider converter
                    rspHlut = matRad_loadHLUT(ct,obj.radiationMode);

                    try
                        % Write Schneider Converter
                        if ~obj.materialConverter.loadConverterFromFile
                            % define density correction
                            matRad_cfg.dispInfo('TOPAS: Writing density correction\n');
                            switch obj.materialConverter.densityCorrection
                                case 'rspHLUT'
                                    densityCorrection.density = [];
                                    for i = 1:size(rspHlut,1)-1
                                        startVal = rspHlut(i,1);
                                        endVal = rspHlut(i+1,1);
                                        range = startVal:1:endVal-1;
                                        densityCorrection.density(end+1:end+numel(range)) = matRad_interp1(rspHlut(:,1),rspHlut(:,2),range);
                                    end
                                    densityCorrection.density(end+1) = rspHlut(end,2); %add last missing value
                                    densityCorrection.boundaries = [rspHlut(1,1) numel(densityCorrection.density)-abs(rspHlut(1,1))];

                                case {'Schneider_TOPAS','Schneider_matRad'}
                                    fname = fullfile(obj.topasFolder,filesep,obj.converterFolder,filesep,obj.infilenames.(['matConv_Schneider_densityCorr_',obj.materialConverter.densityCorrection]));
                                    densityFile = fopen(fname);
                                    densityCorrection.density = fscanf(densityFile,'%f');
                                    fclose(densityFile);
                                    densityCorrection.boundaries = [-1000 numel(densityCorrection.density)-1000];

                            end

                            % define additional density sections
                            switch obj.materialConverter.addSection
                                case 'lung'
                                    addSection = [0.00012 1.05];
                                otherwise
                                    addSection = [];
                            end
                            if exist('addSection','var') && ~isempty(addSection)
                                densityCorrection.density(end+1:end+numel(addSection)) = addSection;
                                densityCorrection.boundaries(end+1) = densityCorrection.boundaries(end)+numel(addSection);
                            end
                            % define Hounsfield Unit Sections
                            switch obj.materialConverter.HUSection
                                case 'default'
                                    densityCorrection.unitSections = [densityCorrection.boundaries];
                                    densityCorrection.offset = 1;
                                    densityCorrection.factor = 0;
                                    densityCorrection.factorOffset = -rspHlut(1,1);

                                case 'advanced'
                                    densityCorrection.offset = [0.00121 1.018 1.03 1.003 1.017 2.201];
                                    densityCorrection.factor = [0.001029700665188 0.000893 0 0.001169 0.000592 0.0005];
                                    densityCorrection.factorOffset = [1000 0 1000 0 0 -2000];

                                    if isfield(obj.materialConverter,'addTitanium') && obj.materialConverter.addTitanium %Titanium independent of set hounsfield unit!
                                        densityCorrection.density(end+1) = 1.00275;
                                        densityCorrection.boundaries(end+1) = densityCorrection.boundaries(end)+1;
                                        densityCorrection.offset(end+1) = 4.54;
                                        densityCorrection.factor(end+1) = 0;
                                        densityCorrection.factorOffset(end+1) = 0;
                                    end

                                    densityCorrection.unitSections = [densityCorrection.boundaries(1) -98 15 23 101 2001 densityCorrection.boundaries(2:end)];
                            end
                            for i = numel(densityCorrection.offset)+1:numel(densityCorrection.unitSections)-1
                                densityCorrection.offset(i) = 1;
                                densityCorrection.factor(i) = 0;
                                densityCorrection.factorOffset(i) = 0;
                            end

                            % write density correction
                            fprintf(fID,'# -- Density correction\n');
                            fprintf(fID,['dv:Ge/Patient/DensityCorrection = %i',repmat(' %.6g',1,numel(densityCorrection.density)),' g/cm3\n'],numel(densityCorrection.density),densityCorrection.density);
                            fprintf(fID,['iv:Ge/Patient/SchneiderHounsfieldUnitSections = %i',repmat(' %g',1,numel(densityCorrection.unitSections)),'\n'],numel(densityCorrection.unitSections),densityCorrection.unitSections);
                            fprintf(fID,['uv:Ge/Patient/SchneiderDensityOffset = %i',repmat(' %g',1,numel(densityCorrection.offset)),'\n'],numel(densityCorrection.offset),densityCorrection.offset);
                            % this is needed for a custom fprintf format which formats integers i to 'i.' and floats without trailing zeros
                            % TODO: check whether this can be removed -> this is potentially not necessary but was done to mimick the original TOPAS Schneider converter file
                            TOPASisFloat = mod(densityCorrection.factor,1)==0;
                            fprintf(fID,['uv:Ge/Patient/SchneiderDensityFactor = %i ',strjoin(cellstr(char('%1.01f '.*TOPASisFloat' + '%1.15g '.*~TOPASisFloat'))),'\n'],numel(densityCorrection.factor),densityCorrection.factor);
                            TOPASisFloat = mod(densityCorrection.factorOffset,1)==0;
                            fprintf(fID,['uv:Ge/Patient/SchneiderDensityFactorOffset = %i ',strjoin(cellstr(char('%1.01f '.*TOPASisFloat' + '%1.15g '.*~TOPASisFloat'))),'\n'],numel(densityCorrection.factorOffset),densityCorrection.factorOffset);
                            %                         fprintf(fID,'uv:Ge/Patient/SchneiderDensityFactor = 8 0.001029700665188 0.000893 0.0 0.001169 0.000592 0.0005 0.0 0.0\n');
                            %                         fprintf(fID,'uv:Ge/Patient/SchneiderDensityFactorOffset = 8 1000. 0. 1000. 0. 0. -2000. 0. 0.0\n\n');

                            % define HU to material sections
                            matRad_cfg.dispInfo('TOPAS: Writing HU to material sections\n');
                            switch obj.materialConverter.HUToMaterial
                                case 'default'
                                    HUToMaterial.sections = rspHlut(2,1);
                                case 'MCsquare'
                                    HUToMaterial.sections = [-1000 -950 -120 -82 -52 -22 8 19 80 120 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500];
                                case 'advanced'
                                    HUToMaterial.sections = [-950 -120 -83 -53 -23 7 18 80 120 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500];
                            end
                            HUToMaterial.sections = [densityCorrection.boundaries(1) HUToMaterial.sections densityCorrection.boundaries(2:end)];
                            % write HU to material sections
                            %                         fprintf(fID,'i:Ge/Patient/MinImagingValue = %d\n',densityCorrection.boundaries(1));
                            fprintf(fID,['iv:Ge/Patient/SchneiderHUToMaterialSections = %i ',repmat('%d ',1,numel(HUToMaterial.sections)),'\n\n'],numel(HUToMaterial.sections),HUToMaterial.sections);
                            % load defined material based on materialConverter.HUToMaterial

                            fname = fullfile(obj.topasFolder,filesep,obj.converterFolder,filesep,obj.infilenames.matConv_Schneider_definedMaterials.(obj.materialConverter.HUToMaterial));
                            materials = strsplit(fileread(fname),'\n')';
                            switch obj.materialConverter.HUToMaterial
                                case 'default'
                                    fprintf(fID,'%s \n',materials{1:end-1});
                                    ExcitationEnergies = str2double(strsplit(materials{end}(strfind(materials{end},'=')+4:end-3)));
                                    if ~isempty(strfind(lower(obj.materialConverter.addSection),lower('lung')))
                                        fprintf(fID,'uv:Ge/Patient/SchneiderMaterialsWeight%i = 5 0.10404040 0.75656566 0.03131313 0.10606061 0.00202020\n',length(materials)-2);
                                        ExcitationEnergies = [ExcitationEnergies' 75.3];
                                    end
                                    fprintf(fID,['dv:Ge/Patient/SchneiderMaterialMeanExcitationEnergy = %i',repmat(' %.6g',1,numel(ExcitationEnergies)),' eV\n'],numel(ExcitationEnergies),ExcitationEnergies);
                                case 'advanced'
                                    fprintf(fID,'\n%s\n',materials{:});
                                case 'MCsquare'
                                    fprintf(fID,'\n%s\n',materials{:});
                            end

                            switch obj.materialConverter.HUToMaterial
                                case 'advanced'
                                    counter = 25;
                                    if isfield(obj.materialConverter,'addTitanium') && obj.materialConverter.addTitanium
                                        fprintf(fID,'uv:Ge/Patient/SchneiderMaterialsWeight%i = 15 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0',counter);
                                        counter = counter + 1;
                                    end
                                    %                                 fprintf(fID,'uv:Ge/Patient/SchneiderMaterialsWeight%i = 15 0.10404040 0.10606061 0.75656566 0.03131313 0.0 0.00202020 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0',counter);
                                    fprintf(fID,'uv:Ge/Patient/SchneiderMaterialsWeight%i = 15 0.101278 0.102310 0.028650 0.757072 0.000730 0.000800 0.002250 0.002660 0.0 0.000090 0.001840 0.001940 0.0 0.000370 0.000010',counter);
                            end
                        else
                            fname = fullfile(obj.topasFolder,filesep,obj.converterFolder,filesep,obj.infilenames.matConv_Schneider_loadFromFile);
                            converter = fileread(fname);
                            fprintf(fID,'\n%s\n',converter);
                        end

                        % write patient environment
                        matRad_cfg.dispInfo('TOPAS: Writing patient environment\n');
                        fprintf(fID,'\n# -- Patient parameters\n');
                        fprintf(fID,'s:Ge/Patient/Parent="World"\n');
                        fprintf(fID,'s:Ge/Patient/Type = "TsImageCube"\n');
                        fprintf(fID,'b:Ge/Patient/DumpImagingValues = "True"\n');
                        fprintf(fID,'s:Ge/Patient/InputDirectory = "./"\n');
                        fprintf(fID,'s:Ge/Patient/InputFile = "%s"\n',dataFile);
                        fprintf(fID,'s:Ge/Patient/ImagingtoMaterialConverter = "Schneider"\n');
                        fprintf(fID,'i:Ge/Patient/NumberOfVoxelsX = %d\n',ct.cubeDim(2));
                        fprintf(fID,'i:Ge/Patient/NumberOfVoxelsY = %d\n',ct.cubeDim(1));
                        fprintf(fID,'iv:Ge/Patient/NumberOfVoxelsZ = 1 %d\n',ct.cubeDim(3));
                        fprintf(fID,'d:Ge/Patient/VoxelSizeX       = %.3f mm\n',ct.resolution.x);
                        fprintf(fID,'d:Ge/Patient/VoxelSizeY       = %.3f mm\n',ct.resolution.y);
                        fprintf(fID,'dv:Ge/Patient/VoxelSizeZ       = 1 %.3f mm\n',ct.resolution.z);
                        fprintf(fID,'s:Ge/Patient/DataType  = "SHORT"\n');

                        fclose(fID);

                        % write HU data
                        matRad_cfg.dispInfo('TOPAS: Export patient cube\n');
                        huCube = int32(permute(ct.cubeHU{ctScen},permutation));
                        fID = fopen(fullfile(obj.workingDir, dataFile),'w');
                        fwrite(fID,huCube,'short');
                        fclose(fID);
                        cube = huCube;
                    catch ME
                        matRad_cfg.dispWarning(ME.message);
                        matRad_cfg.dispError(['TOPAS: Error in Schneider Converter! (line ',num2str(ME.stack(1).line),')']);
                    end

                otherwise
                    matRad_cfg.dispError('Material Conversion rule "%s" not implemented (yet)!\n',obj.materialConverter.mode);
            end
            obj.MCparam.imageCube{ctScen} = cube;


        end

        function writeRangeShifter(~,fID,rangeShifter,sourceToNozzleDistance)
            %TODO: Insert documentation
            %Hardcoded PMMA range shifter for now
            pmma_rsp = 1.165;
            rsWidth = rangeShifter.eqThickness / pmma_rsp;

            fprintf(fID,'s:Ge/%s/Parent   = "Nozzle"\n',rangeShifter.topasID);
            fprintf(fID,'s:Ge/%s/Type     = "TsBox"\n',rangeShifter.topasID);
            fprintf(fID,'s:Ge/%s/Material = "Lucite"\n',rangeShifter.topasID);
            fprintf(fID,'d:Ge/%s/HLX      = 250 mm\n',rangeShifter.topasID);
            fprintf(fID,'d:Ge/%s/HLY      = 250 mm\n',rangeShifter.topasID);
            fprintf(fID,'d:Ge/%s/HLZ      = %f  mm\n',rangeShifter.topasID,rsWidth/2);
            fprintf(fID,'d:Ge/%s/TransX   = 500 mm * Tf/Beam/%sOut/Value\n',rangeShifter.topasID,rangeShifter.topasID);
            fprintf(fID,'d:Ge/%s/TransY   = 0   mm\n',rangeShifter.topasID);
            fprintf(fID,'d:Ge/%s/TransZ   = %f mm\n',rangeShifter.topasID,rangeShifter.sourceRashiDistance - sourceToNozzleDistance);

        end

        function writeMCparam(obj)
            %TODO: Insert documentation
            %write MCparam file with basic parameters
            MCparam = obj.MCparam;
            save(fullfile(obj.workingDir,'MCparam.mat'),'MCparam','-v7');
        end

    end
    methods(Static)
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
                checkModality = any(strcmp(DoseEngines.matRad_TopasMCEngine.possibleRadiationModes, machine.meta.radiationMode));
                
                preCheck = checkBasic && checkModality;

                if ~preCheck
                    return;
                else
                    available = preCheck;
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return;
            end

            
       end
    end
end

