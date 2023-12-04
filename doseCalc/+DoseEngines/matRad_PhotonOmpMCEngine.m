classdef matRad_PhotonOmpMCEngine < DoseEngines.matRad_MonteCarloEngineAbstract
    % Engine for photon dose calculation based on monte carlo
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
        possibleRadiationModes = 'photons';
        name = 'ompMC';
        shortName = 'ompMC';
    end

    properties (SetAccess = public, GetAccess = public)
        visBool = false; %binary switch to en/disable visualitzation
        useCornersSCD = true; %false -> use ISO corners
        
         % This factor calibrates to 1 Gy in a %(5x5)cm^2 open field (1 bixel) at
         % 5cm depth for SSD = 900 which corresponds to the calibration for the
         % analytical base data.
        absCalibrationFactor = 3.49056 * 1e12; %Approximate!
        
        omcFolder;
    end

    properties (SetAccess = protected, GetAccess = public)
        ompMCoptions;
        ompMCgeo;
        ompMCsource;

        cubeHU;         %resamples HU cube
        cubeRho;        %density cube
        cubeMatIx;      %material assignment
    end

    properties (Constant)
        scale = 10;
    end

    methods
        function this = matRad_PhotonOmpMCEngine(pln)
            % Constructor
            %
            % call
            %   engine = DoseEngines.matRad_DoseEnginePhotonsOmpMCct,stf,pln,cst)
            %
            % input
            %   pln:                        matRad plan meta information struct

            if nargin < 1
                pln = [];
            end

            % call superclass constructor
            this = this@DoseEngines.matRad_MonteCarloEngineAbstract(pln);

            matRad_cfg = MatRad_Config.instance();
            this.omcFolder = [matRad_cfg.matRadRoot filesep 'ompMC'];

            if ~matRad_checkMexFileExists('omc_matrad') %exist('matRad_ompInterface','file') ~= 3
                matRad_cfg.dispWarning('Compiled mex interface not found. Trying to compile the ompMC interface on the fly!');
                try
                    this.compileOmpMCInterface(this.omcFolder);
                catch MException
                    matRad_cfg.dispError('Could not find/generate mex interface for MC dose calculation.\nCause of error:\n%s\n Please compile it yourself (preferably with OpenMP support).',MException.message);
                end
            end
        end
    end

    methods (Access = protected)
        function dij = calcDose(this,ct,cst,stf)
            % matRad ompMC monte carlo photon dose calculation wrapper
            % can be automaticly called through matRad_calcDose or
            % matRad_calcPhotonDoseMC
            %
            % call
            %   dij = this.calcDose(ct,stf,pln,cst)
            %
            % input
            %   ct:                         matRad ct struct
            %   stf:                        matRad steering information struct
            %   cst:                        matRad cst struct
            % output
            %   dij:                        matRad dij struct
            %
            % References
            %   -
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Copyright 2018 the matRad development team.
            %
            % This file is part of the matRad project. It is subject to the license
            % terms in the LICENSE file found in the top-level directory of this
            % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
            % of the matRad project, including this file, may be copied, modified,
            % propagated, or distributed except according to the terms contained in the
            % LICENSE file.
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            matRad_cfg =  MatRad_Config.instance();

            %run initDoseCalc as usual
            dij = this.initDoseCalc(ct,cst,stf);

            %ompMC for matRad returns dose/history * nHistories.
               
            bixelWidth = unique([stf.bixelWidth]);

            if numel(bixelWidth) > 1
                matRad_cfg.dispWarning('Varying bixel width in stf, calibartion might be wrong!')
                bixelWidth = mean(bixelWidth);
            end           

            %Now we have to calibrate to the the beamlet width.
            calibrationFactor = this.absCalibrationFactor * (bixelWidth/50)^2;

            scenCount = 0;
            %run over all scenarios
            for scenarioIx = 1:this.multScen.totNumScen
                ctScen = this.multScen.linearMask(scenarioIx,1);
                shiftScen = this.multScen.linearMask(scenarioIx,2);
                rangeShiftScen = this.multScen.linearMask(scenarioIx,3);

                if this.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                    scenCount = scenCount + 1;
                    
                    % manipulate isocenter
                    shiftedIsoCenter = vertcat(stf(:).isoCenter) + this.multScen.isoShift(scenarioIx,:) + dij.doseGrid.isoCenterOffset;
                    this.ompMCgeo.isoCenter = shiftedIsoCenter;
                    tmpStf = stf;
                    
                    for k = 1:length(tmpStf)
                        tmpStf(k).isoCenter = shiftedIsoCenter;
                    end

                    % load ompMC source
                    this.getOmpMCsource(tmpStf);

                    % Book keeping for dij
                    counter = 0;
                    for i = 1:dij.numOfBeams
                        for j = 1:tmpStf(i).numOfRays
                            counter = counter + 1;
                            dij.beamNum(counter)  = i;
                            dij.rayNum(counter)   = j;
                            dij.bixelNum(counter) = j;
                        end
                    end

                    if this.multScen.totNumScen == 1
                        matRad_cfg.dispInfo('matRad: OmpMC photon dose calculation... \n');
                    else
                        matRad_cfg.dispInfo('matRad: OmpMC photon dose calculation for scenario %d of %d... \n',scenCount,this.multScen.totNumScen);
                    end

                    %Call the Monte Carlo simulation and catch  possible mex
                    %interface issues
                    try
                        %If we ask for variance, a field in the dij will be filled
                        if this.outputMCvariance
                            [dij.physicalDose{ctScen,shiftScen,rangeShiftScen},dij.physicalDose_MCvar{ctScen,shiftScen,rangeShiftScen}] = omc_matrad(this.cubeRho{ctScen},this.cubeMatIx{ctScen},this.ompMCgeo,this.ompMCsource,this.ompMCoptions);
                        else
                            [dij.physicalDose{ctScen,shiftScen,rangeShiftScen}] = omc_matrad(this.cubeRho{ctScen},this.cubeMatIx{ctScen},this.ompMCgeo,this.ompMCsource,this.ompMCoptions);
                        end
                    catch ME
                        errorString = [ME.message '\nThis error was thrown by the MEX-interface of ompMC.\nMex interfaces can raise compatability issues which may be resolved by compiling them by hand directly on your particular system.'];
                        matRad_cfg.dispError(errorString);
                    end

                    %Calibrate the dose with above factor
                    dij.physicalDose{scenarioIx} = dij.physicalDose{scenarioIx} * calibrationFactor;
                    if isfield(dij,'physicalDose_MCvar')
                        dij.physicalDose_MCvar{scenarioIx} = dij.physicalDose_MCvar{scenarioIx} * calibrationFactor^2;
                    end
                end
            end

            %Finalize dose calculation
            dij = this.finalizeDose(dij);   
        end
    
        function dij = initDoseCalc(this,ct,cst,stf)

            dij = initDoseCalc@DoseEngines.matRad_MonteCarloEngineAbstract(this,ct,cst,stf);
            
            matRad_cfg = MatRad_Config.instance();

            % set up arrays for book keeping
            dij.bixelNum = NaN*ones(dij.totalNumOfBixels,1);
            dij.rayNum   = NaN*ones(dij.totalNumOfBixels,1);
            dij.beamNum  = NaN*ones(dij.totalNumOfBixels,1);

            dij.numHistoriesPerBeamlet = this.numHistoriesPerBeamlet;

            %% Setup OmpMC options / parameters
            this.setOmpMCoptions();

            % conversion from HU to densities & materials
            this.materialConversion(dij.ctGrid,dij.doseGrid,ct);
            
            % Create the Geometry
            this.getOmpMCgeometry(dij.doseGrid);

            %% Create beamlet source
            this.getOmpMCsource(stf);
        end
    end

    methods (Access = private)
        function setOmpMCoptions(obj)
            matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class
            
            %display options
            obj.ompMCoptions.verbose = matRad_cfg.logLevel - 1;
            
            % start MC control
            obj.ompMCoptions.nHistories = obj.numHistoriesPerBeamlet;
            obj.ompMCoptions.nSplit = 20;
            obj.ompMCoptions.nBatches = 10;
            obj.ompMCoptions.randomSeeds = [97 33];
            
            %start source definition
            obj.ompMCoptions.spectrumFile       = [obj.omcFolder filesep 'spectra' filesep 'mohan6.spectrum'];
            obj.ompMCoptions.monoEnergy         = 0.1;
            obj.ompMCoptions.charge             = 0;
            obj.ompMCoptions.sourceGeometry     = 'gaussian';
            obj.ompMCoptions.sourceGaussianWidth = obj.getSourceWidthFromPenumbra./obj.scale;
            
            % start MC transport
            obj.ompMCoptions.dataFolder   = [obj.omcFolder filesep 'data' filesep];
            obj.ompMCoptions.pegsFile     = [obj.omcFolder filesep 'pegs4' filesep '700icru.pegs4dat'];
            obj.ompMCoptions.pgs4formFile = [obj.omcFolder filesep 'pegs4' filesep 'pgs4form.dat'];
            
            obj.ompMCoptions.global_ecut = 0.7;
            obj.ompMCoptions.global_pcut = 0.010;
            
            % Relative Threshold for dose
            obj.ompMCoptions.relDoseThreshold = 1 - obj.relativeDosimetricCutOff;
            
            % Output folders
            obj.ompMCoptions.outputFolder = [obj.omcFolder filesep 'output' filesep];            
        end

        function sigmaGauss = getSourceWidthFromPenumbra(obj)
            % gaussian filter to model penumbra from (measured) machine output / see diploma thesis siggel 4.1.2
            matRad_cfg = MatRad_Config.instance();
            
            if isfield(obj.machine.data,'penumbraFWHMatIso')
                penumbraFWHM = obj.machine.data.penumbraFWHMatIso;
            else
                penumbraFWHM = 5;
                matRad_cfg.dispWarning('photon machine file does not contain measured penumbra width in machine.data.penumbraFWHMatIso. Assuming 5 mm.');
            end
            
            sourceFWHM = penumbraFWHM * obj.machine.meta.SCD/(obj.machine.meta.SAD - obj.machine.meta.SCD);
            sigmaGauss = sourceFWHM / sqrt(8*log(2)); % [mm]  
        end

        function obj = materialConversion(obj,ctGrid,doseGrid,ct)
            % conversion from HU to densities & materials
            obj.cubeHU      = cell(1,ct.numOfCtScen);
            obj.cubeMatIx   = cell(1,ct.numOfCtScen);
            obj.cubeRho     = cell(1,ct.numOfCtScen);
            
            % Create Material Density Cube
            material = obj.setupMaterials();
            
            for s = 1:ct.numOfCtScen
                
                obj.cubeHU{s} =  matRad_interp3(ctGrid.x,ctGrid.y',ctGrid.z,ct.cubeHU{s}, ...
                    doseGrid.x,doseGrid.y',doseGrid.z,'nearest');
                
                % projecting out of bounds HU values where necessary
                if max(obj.cubeHU{s}(:)) > material{end,3}
                    matRad_cfg.dispWarning('Projecting out of range HU values');
                    obj.cubeHU{s}(obj.cubeHU{s}(:) > material{end,3}) = material{end,3};
                end
                if min(obj.cubeHU{s}(:)) < material{1,2}
                    matRad_cfg.dispWarning('Projecting out of range HU values');
                    obj.cubeHU{s}(obj.cubeHU{s}(:) < material{1,2}) = material{1,2};
                end
                
                % find material index
                obj.cubeMatIx{s} = NaN*ones(doseGrid.dimensions,'int32');
                for i = size(material,1):-1:1
                    obj.cubeMatIx{s}(obj.cubeHU{s} <= material{i,3}) = i;
                end
                
                % create an artificial HU lookup table
                hlut = [];
                for i = 1:size(material,1)
                    hlut = [hlut;material{i,2} material{i,4};material{i,3}-1e-10 material{i,5}]; % add eps for interpolation
                end
                
                obj.cubeRho{s} = interp1(hlut(:,1),hlut(:,2),obj.cubeHU{s});
                
            end
        end

        function getOmpMCgeometry(obj,doseGrid)
            obj.ompMCgeo.xBounds = (doseGrid.resolution.y * (0.5 + [0:doseGrid.dimensions(1)])) ./ obj.scale;
            obj.ompMCgeo.yBounds = (doseGrid.resolution.x * (0.5 + [0:doseGrid.dimensions(2)])) ./ obj.scale;
            obj.ompMCgeo.zBounds = (doseGrid.resolution.z * (0.5 + [0:doseGrid.dimensions(3)])) ./ obj.scale;
        
            % Create Material Density Cube
            obj.ompMCgeo.material    = obj.setupMaterials();

            %% debug visualization
            if obj.visBool

                figure;
                hold on;

                axis equal;

                % ct box
                ctCorner1 = [obj.ompMCgeo.xBounds(1) obj.ompMCgeo.yBounds(1) obj.ompMCgeo.zBounds(1)];
                ctCorner2 = [obj.ompMCgeo.xBounds(end) obj.ompMCgeo.yBounds(end) obj.ompMCgeo.zBounds(end)];
                plot3([ctCorner1(1) ctCorner2(1)],[ctCorner1(2) ctCorner1(2)],[ctCorner1(3) ctCorner1(3)],'k' );
                plot3([ctCorner1(1) ctCorner2(1)],[ctCorner2(2) ctCorner2(2)],[ctCorner1(3) ctCorner1(3)],'k' );
                plot3([ctCorner1(1) ctCorner1(1)],[ctCorner1(2) ctCorner2(2)],[ctCorner1(3) ctCorner1(3)],'k' );
                plot3([ctCorner2(1) ctCorner2(1)],[ctCorner1(2) ctCorner2(2)],[ctCorner1(3) ctCorner1(3)],'k' );
                plot3([ctCorner1(1) ctCorner2(1)],[ctCorner1(2) ctCorner1(2)],[ctCorner2(3) ctCorner2(3)],'k' );
                plot3([ctCorner1(1) ctCorner2(1)],[ctCorner2(2) ctCorner2(2)],[ctCorner2(3) ctCorner2(3)],'k' );
                plot3([ctCorner1(1) ctCorner1(1)],[ctCorner1(2) ctCorner2(2)],[ctCorner2(3) ctCorner2(3)],'k' );
                plot3([ctCorner2(1) ctCorner2(1)],[ctCorner1(2) ctCorner2(2)],[ctCorner2(3) ctCorner2(3)],'k' );
                plot3([ctCorner1(1) ctCorner1(1)],[ctCorner1(2) ctCorner1(2)],[ctCorner1(3) ctCorner2(3)],'k' );
                plot3([ctCorner2(1) ctCorner2(1)],[ctCorner1(2) ctCorner1(2)],[ctCorner1(3) ctCorner2(3)],'k' );
                plot3([ctCorner1(1) ctCorner1(1)],[ctCorner2(2) ctCorner2(2)],[ctCorner1(3) ctCorner2(3)],'k' );
                plot3([ctCorner2(1) ctCorner2(1)],[ctCorner2(2) ctCorner2(2)],[ctCorner1(3) ctCorner2(3)],'k' );

                xlabel('x [cm]');
                ylabel('y [cm]');
                zlabel('z [cm]');

                rotate3d on;

            end
        end

        function material = setupMaterials(~)
            material = cell(4,5);
            material{1,1} = 'AIR700ICRU';
            material{1,2} = -1024;
            material{1,3} = -974;
            material{1,4} = 0.001;
            material{1,5} = 0.044;
            material{2,1} = 'LUNG700ICRU';
            material{2,2} = -974;
            material{2,3} = -724;
            material{2,4} = 0.044;
            material{2,5} = 0.302;
            material{3,1} = 'ICRUTISSUE700ICRU';
            material{3,2} = -724;
            material{3,3} = 101;
            material{3,4} = 0.302;
            material{3,5} = 1.101;
            material{4,1} = 'ICRPBONE700ICRU';
            material{4,2} = 101;
            material{4,3} = 1976;
            material{4,4} = 1.101;
            material{4,5} = 2.088;
            
        end
        
        function getOmpMCsource(obj,stf)
            numOfBeams = numel(stf);

            numOfBixelsPerBeam = [stf(:).numOfRays];
            totalNumOfBixels = sum(numOfBixelsPerBeam);
            beamSource = zeros(numOfBeams, 3);

            bixelCorner = zeros(totalNumOfBixels,3);
            bixelSide1 = zeros(totalNumOfBixels,3);
            bixelSide2 = zeros(totalNumOfBixels,3);

            beamNum = zeros(1,prod(totalNumOfBixels,numOfBeams));
            counter = 0;

            for i = 1:numOfBeams % loop over all beams

                % define beam source in physical coordinate system in cm
                beamSource(i,:) = (stf(i).sourcePoint + stf(i).isoCenter)/10;

                for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!

                    counter = counter + 1;

                    beamNum(counter)  = i;

                    if obj.useCornersSCD
                        beamletCorners = stf(i).ray(j).rayCorners_SCD;
                    else
                        beamletCorners = stf(i).ray(j).beamletCornersAtIso;
                    end

                    % get bixel corner and delimiting vectors.
                    % a) change coordinate system (Isocenter cs-> physical cs) and units mm -> cm
                    currCorner = (beamletCorners(1,:) + stf(i).isoCenter) ./ obj.scale;
                    bixelCorner(counter,:) = currCorner;
                    bixelSide1(counter,:) = (beamletCorners(2,:) + stf(i).isoCenter) ./ obj.scale - currCorner;
                    bixelSide2(counter,:) = (beamletCorners(4,:) + stf(i).isoCenter) ./ obj.scale - currCorner;

                    if obj.visBool
                        for k = 1:4
                            currCornerVis = (beamletCorners(k,:) + stf(i).isoCenter)/10;
                            % rays connecting source and ray corner
                            plot3([beamSource(i,1) currCornerVis(1)],[beamSource(i,2) currCornerVis(2)],[beamSource(i,3) currCornerVis(3)],'b')
                            % connection between corners
                            lRayCorner = (beamletCorners(mod(k,4) + 1,:) + stf(i).isoCenter)/10;
                            plot3([lRayCorner(1) currCornerVis(1)],[lRayCorner(2) currCornerVis(2)],[lRayCorner(3) currCornerVis(3)],'r')
                        end
                    end

                end

            end

            obj.ompMCsource.nBeams = numOfBeams;
            obj.ompMCsource.iBeam = beamNum(:);

            % Switch x and y directions to match ompMC cs.
            obj.ompMCsource.xSource = beamSource(:,2);
            obj.ompMCsource.ySource = beamSource(:,1);
            obj.ompMCsource.zSource = beamSource(:,3);

            obj.ompMCsource.nBixels = sum(numOfBixelsPerBeam(:));
            obj.ompMCsource.xCorner = bixelCorner(:,2);
            obj.ompMCsource.yCorner = bixelCorner(:,1);
            obj.ompMCsource.zCorner = bixelCorner(:,3);

            obj.ompMCsource.xSide1 = bixelSide1(:,2);
            obj.ompMCsource.ySide1 = bixelSide1(:,1);
            obj.ompMCsource.zSide1 = bixelSide1(:,3);

            obj.ompMCsource.xSide2 = bixelSide2(:,2);
            obj.ompMCsource.ySide2 = bixelSide2(:,1);
            obj.ompMCsource.zSide2 = bixelSide2(:,3);

            if obj.visBool
                plot3(obj.ompMCsource.ySource,obj.ompMCsource.xSource,obj.ompMCsource.zSource,'rx')
            end

        end
    
    end

    methods (Static)

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
                checkModality = any(strcmp(DoseEngines.matRad_PhotonOmpMCEngine.possibleRadiationModes, machine.meta.radiationMode));

                preCheck = checkBasic && checkModality;

                if ~preCheck
                    return;
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return;
            end

            checkMeta = all(isfield(machine.meta,{'SAD','SCD'}));

            if checkMeta
                available = true;
                msg = 'The ompMC machine is not representing the machine exactly and approximates it with a virtual Gaussian source and generic primary fluence & 6 MV energy spectrum!';
            end
        end

        function compileOmpMCInterface(dest,omcFolder)
            % Compiles the ompMC interface (integrated as submodule)
            %
            % call
            %   matRad_OmpConfig.compileOmpMCInterface()
            %   matRad_OmpConfig.compileOmpMCInterface(dest)
            %   matRad_OmpConfig.compileOmpMCInterface(dest,sourceFolder)
            % if an object is instantiated, matRad_OmpConfig can be replaced by the
            % object handle
            %
            % input:
            %   dest:           (optional) destination for mex file. Default: location
            %                   of this file
            %   sourceFolder:   (optional) path to ompMC . Default assumes its checked
            %                   out in the submodules folder of matRad
            %
            % References
            %
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
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            matRad_cfg = MatRad_Config.instance();
            
            env = matRad_getEnvironment();
            
            %Our destination usually lies in the ompMC thirdPartyFolder
            if nargin < 1
                dest = [matRad_cfg.matRadRoot filesep 'ompMC'];
            end
            
            %We can recompile form the submodules
            if nargin < 2
                omcFolder = [matRad_cfg.matRadRoot filesep 'submodules' filesep 'ompMC'];
            end
            
            sourceFolder = [omcFolder filesep 'src'];
            interfaceFolder = [omcFolder filesep 'ucodes' filesep 'omc_matrad'];
            
            mainFile = [interfaceFolder filesep 'omc_matrad.c'];
            
            addFiles = {'ompmc.c','omc_utilities.c','omc_random.c'};
            addFiles = cellfun(@(f) fullfile(sourceFolder,f),addFiles,'UniformOutput',false);
            
            addFiles = strjoin(addFiles,' ');
            
            if exist ('OCTAVE_VERSION','builtin')
                ccName = evalc('mkoctfile -p CC');
            else
                myCCompiler = mex.getCompilerConfigurations('C','Selected');
                ccName = myCCompiler.ShortName;
            end
            
            %These settings have only been tested for MSVC and g++. You may need to adapt for other compilers
            if ~isempty(strfind(ccName,'MSVC')) %Not use contains(...) because of octave
                flags{1,1} = 'COMPFLAGS';
                flags{1,2} = '/openmp';
                flags{2,1} = 'OPTIMFLAGS';
                flags{2,2} = '/O2';
            else
                flags{1,1} = 'CFLAGS';
                flags{1,2} = '-std=gnu99 -fopenmp -O3';
                flags{2,1} = 'LDFLAGS';
                flags{2,2} = '-fopenmp';
                
            end
            
            includestring =  ['-I' sourceFolder];
            
            flagstring = '';
            
            %For Octave, the flags will be set in the environment, while they
            %will be parsed as string arguments in MATLAB
            for flag = 1:size(flags,1)
                if strcmp(env,'OCTAVE')
                    preFlagContent = eval(['mkoctfile -p ' flags{flag,1}]);
                    if ~isempty(preFlagContent)
                        preFlagContent = preFlagContent(1:end-1); %Strip newline
                    end
                    newContent = [preFlagContent ' ' flags{flag,2}];
                    setenv(flags{flag,1},newContent);
                    matRad_cfg.dispDebug('Set compiler flag %s to %s\n',flags{flag,1},newContent);
                else
                    flagstring = [flagstring flags{flag,1} '="$' flags{flag,1} ' ' flags{flag,2} '" '];
                end
            end
            
            mexCall = ['mex -largeArrayDims ' flagstring ' ' includestring ' ' mainFile ' ' addFiles];
            matRad_cfg.dispDebug('Compiler call: %s\n',mexCall);
            
            currDir = pwd;
            cd(dest);
            eval(mexCall);
            cd(currDir);
        end
    end
end

