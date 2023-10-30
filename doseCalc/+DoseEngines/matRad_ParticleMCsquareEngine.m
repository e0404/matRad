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
            
            % call superclass constructor
            this = this@DoseEngines.matRad_MonteCarloEngineAbstract(pln);

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


            %% check if binaries are available
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
            fullfilename = mfilename('fullpath');
            MCsquareFolder = [matRad_cfg.matRadRoot filesep 'MCsquare' filesep 'bin'];

            % cd to MCsquare folder (necessary for binary)
            cd(MCsquareFolder);

            %Check Materials
            if ~exist([MCsquareFolder filesep 'Materials'],'dir') || ~exist(fullfile(MCsquareFolder,'Materials','list.dat'),'file')
                matRad_cfg.dispInfo('First call of MCsquare: unzipping Materials...');    
                unzip('Materials.zip');
                matRad_cfg.dispInfo('Done');
            end


            %Now we can run calcDoseInit as usual
            [dij,ct,cst,stf] = this.calcDoseInit(ct,cst,stf);
            
            % We need to adjust the offset used in matRad_calcDoseInit
            mcSquareAddIsoCenterOffset = [dij.doseGrid.resolution.x/2 dij.doseGrid.resolution.y/2 dij.doseGrid.resolution.z/2] ...
                            - [dij.ctGrid.resolution.x   dij.ctGrid.resolution.y   dij.ctGrid.resolution.z];
            mcSquareAddIsoCenterOffset = mcSquareAddIsoCenterOffset - dij.doseGrid.isoCenterOffset;

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

            MCsquareConfig = MatRad_MCsquareConfig;

            MCsquareConfig.BDL_Plan_File = 'currBixels.txt';
            MCsquareConfig.CT_File       = 'MC2patientCT.mhd';
            MCsquareConfig.Num_Threads   = this.nbThreads;
            MCsquareConfig.RNG_Seed      = 1234;
            MCsquareConfig.Num_Primaries = this.numHistoriesPerBeamlet;

            % turn simulation of individual beamlets
            MCsquareConfig.Beamlet_Mode = ~this.calcDoseDirect;
            % turn of writing of full dose cube
            MCsquareConfig.Dose_MHD_Output = this.calcDoseDirect;
            % turn on sparse output
            MCsquareConfig.Dose_Sparse_Output = ~this.calcDoseDirect;
            % set threshold of sparse matrix generation
            MCsquareConfig.Dose_Sparse_Threshold = 1 - this.relativeDosimetricCutOff;

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
                            if this.calcDoseDirect
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
            [status,cmdout] = system([this.mcSquareBinary ' ' MCsquareConfigFile],'-echo');

            mask = false(dij.doseGrid.numOfVoxels,1);
            mask(this.VdoseGrid) = true;

            % read sparse matrix
            if ~this.calcDoseDirect
                dij.physicalDose{1} = absCalibrationFactorMC2 * matRad_sparseBeamletsReaderMCsquare ( ...
                                [MCsquareConfig.Output_Directory filesep 'Sparse_Dose.bin'], ...
                                dij.doseGrid.dimensions, ...
                                dij.totalNumOfBixels, ...
                                mask);
            else
                cube = matRad_readMhd(MCsquareConfig.Output_Directory,'Dose.mhd');
                dij.physicalDose{1} = sparse(this.VdoseGrid,ones(numel(this.VdoseGrid),1), ...
                                             absCalibrationFactorMC2 * cube(this.VdoseGrid), ...
                                             dij.doseGrid.numOfVoxels,1);
            end

            % reorder influence matrix to comply with matRad default ordering
            if MCsquareConfig.Beamlet_Mode
                dij.physicalDose{1} = dij.physicalDose{1}(:,MCsquareOrder);            
            end        

            matRad_cfg.dispInfo('matRad: done!\n');

            %% clear all data
            %could also be moved to the "finalize" function
            delete([MCsquareConfig.CT_File(1:end-4) '.*']);
            delete('currBixels.txt');
            delete('MCsquareConfig.txt');

            %For Octave temporarily disable confirmation for recursive rmdir
            if strcmp(matRad_cfg.env,'OCTAVE')    
                rmdirConfirmState = confirm_recursive_rmdir(0);
            end
            rmdir(MCsquareConfig.Output_Directory,'s');

            %Reset to old confirmatoin state
            if strcmp(matRad_cfg.env,'OCTAVE')
                confirm_recursive_rmdir(rmdirConfirmState);
            end

            % cd back
            cd(this.currFolder);

            %Finalize dose calculation
            dij = this.calcDoseFinalize(ct,cst,stf,dij);

        end
        
    end
    
    methods(Access = protected)
              
        function setUp(this,nCasePerBixel,calcDoseDirect)    
        % SETUP Set up properties used for dose calculation
        %
        % input:
        %   nCasePerBixel:  number of histories per beamlet
        %   calcDoseDirect: binary switch to enable forward dose calculation output
        %
   
            matRad_cfg = MatRad_Config.instance();

            % first argument should be nCasePerBixel
            if (exist('nCasePerBixel','var') && isa(nCasePerBixel,'numeric'))    
                this.nCasePerBixel = nCasePerBixel;
            else
                %set number of particles simulated per pencil beam
                this.nCasePerBixel = matRad_cfg.propMC.MCsquare_defaultHistories;
                matRad_cfg.dispInfo('No number of Histories given or wrong type given. Using default number of Histories per Bixel: %d\n',this.nCasePerBixel);
            end

            if (exist('calcDoseDirect', 'var'))
                this.calcDoseDirect = true;
            end       
            
        end
        
        function setBinaries(this)
            % setBinaries check if the binaries are available on the current
            % machine and sets to the mcsquarebinary object property
            %

            [~,binaryFile] = this.checkBinaries();
            this.mcSquareBinary = binaryFile;
        end
        
        function [dij,ct,cst,stf] = calcDoseInit(this,ct,cst,stf)
            %% Assingn and check parameters
            
            matRad_cfg = MatRad_Config.instance();            

            %% Call Superclass init function
            [dij,ct,cst,stf] = calcDoseInit@DoseEngines.matRad_MonteCarloEngineAbstract(this,ct,cst,stf); 

            % Since MCsquare 1.1 only allows similar resolution in x&y, we do some
            % extra checks on that before calling the normal calcDoseInit. First, we make sure a
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
            
            if ~hasBinaries
                return;
            end

            %MCsquare currently only works for the generic_MCSquare machine
            available = strcmp(pln.machine,'generic_MCsquare');
            msg = 'Machine check is currently not reliable, machine only works reliably with generic_MCsquare';
        end
              
    end
end

