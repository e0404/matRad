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
    end

    properties (SetAccess = public, GetAccess = public)
        visBool = false; %binary switch to en/disable visualitzation
        useCornersSCD = true; %false -> use ISO corners

        omcFolder;
    end

    properties (SetAccess = protected, GetAccess = public)
        ompMCoptions;
        ompMCgeo;
        ompMCsource;
        cubeRho;        %density cube
        cubeMatIx;      %material assignment
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

            % call superclass constructor
            this = this@DoseEngines.matRad_MonteCarloEngineAbstract(pln);

            matRad_cfg = MatRad_Config.instance();
            this.omcFolder = [matRad_cfg.matRadRoot filesep 'ompMC'];

            fileFolder = fileparts(mfilename('fullpath'));
            if ~matRad_checkMexFileExists('omc_matrad') %exist('matRad_ompInterface','file') ~= 3
                matRad_cfg.dispWarning('Compiled mex interface not found. Trying to compile the ompMC interface on the fly!');
                try
                    matRad_compileOmpMCInterface();
                catch MException
                    matRad_cfg.dispError('Could not find/generate mex interface for MC dose calculation.\nCause of error:\n%s\n Please compile it yourself (preferably with OpenMP support).',MException.message);
                end
            end
        end

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

            %run calcDoseInit as usual
            [dij,ct,cst,stf] = this.calcDoseInit(ct,cst,stf);

            %% Call the OmpMC interface

            %ompMC for matRad returns dose/history * nHistories.
            % This factor calibrates to 1 Gy in a %(5x5)cm^2 open field (1 bixel) at
            % 5cm depth for SSD = 900 which corresponds to the calibration for the
            % analytical base data.
            absCalibrationFactor = 3.49056 * 1e12; %Approximate!
            
            bixelWidth = unique([stf.bixelWidth]);

            if numel(bixelWidth) > 1
                matRad_cfg.dispWarning('Varying bixel width in stf, calibartion might be wrong!')
                bixelWidth = mean(bixelWidth);
            end

            %Now we have to calibrate to the the beamlet width.
            absCalibrationFactor = absCalibrationFactor * (bixelWidth/50)^2;

            %run over all scenarios
            for s = 1:dij.numOfScenarios
                this.ompMCgeo.isoCenter = [stf(:).isoCenter];

                %Run the Monte Carlo simulation and catch  possible mex-interface
                %issues
                try
                    %If we ask for variance, a field in the dij will be filled
                    if this.outputMCvariance
                        [dij.physicalDose{s},dij.physicalDose_MCvar{s}] = omc_matrad(this.cubeRho{s},this.cubeMatIx{s},this.ompMCgeo,this.ompMCsource,this.ompMCoptions);
                    else
                        [dij.physicalDose{s}] = omc_matrad(this.cubeRho{s},this.cubeMatIx{s},this.ompMCgeo,this.ompMCsource,this.ompMCoptions);
                    end
                catch ME
                    errorString = [ME.message '\nThis error was thrown by the MEX-interface of ompMC.\nMex interfaces can raise compatability issues which may be resolved by compiling them by hand directly on your particular system.'];
                    matRad_cfg.dispError(errorString);
                end

                %Calibrate the dose with above factor
                dij.physicalDose{s} = dij.physicalDose{s} * absCalibrationFactor;
                if isfield(dij,'physicalDose_MCvar')
                    dij.physicalDose_MCvar{s} = dij.physicalDose_MCvar{s} * absCalibrationFactor^2;
                end
            end

            %Finalize dose calculation
            dij = this.calcDoseFinalize(ct,cst,stf,dij);   
        end
    end

    methods (Access = protected)
        function [dij,ct,cst,stf] = calcDoseInit(this,ct,cst,stf)

            [dij,ct,cst,stf] = calcDoseInit@DoseEngines.matRad_MonteCarloEngineAbstract(this,ct,cst,stf);
            
            matRad_cfg = MatRad_Config.instance();

            % gaussian filter to model penumbra from (measured) machine output / see diploma thesis siggel 4.1.2
            if isfield(this.machine.data,'penumbraFWHMatIso')
                penumbraFWHM = this.machine.data.penumbraFWHMatIso;
            else
                penumbraFWHM = 5;
                matRad_cfg.dispWarning('photon machine file does not contain measured penumbra width in machine.data.penumbraFWHMatIso. Assuming 5 mm.');
            end

            sourceFWHM = penumbraFWHM * this.machine.meta.SCD/(this.machine.meta.SAD - this.machine.meta.SCD);
            sigmaGauss = sourceFWHM / sqrt(8*log(2)); % [mm]

            % set up arrays for book keeping
            dij.bixelNum = NaN*ones(dij.totalNumOfBixels,1);
            dij.rayNum   = NaN*ones(dij.totalNumOfBixels,1);
            dij.beamNum  = NaN*ones(dij.totalNumOfBixels,1);

            dij.numHistoriesPerBeamlet = this.numHistoriesPerBeamlet;

            %% Setup OmpMC options / parameters

            %display options
            this.ompMCoptions.verbose = matRad_cfg.logLevel - 1;

            % start MC control
            this.ompMCoptions.nHistories = this.numHistoriesPerBeamlet;
            this.ompMCoptions.nSplit = 20;
            this.ompMCoptions.nBatches = 10;
            this.ompMCoptions.randomSeeds = [97 33];

            %start source definition
            this.ompMCoptions.spectrumFile       = [this.omcFolder filesep 'spectra' filesep 'mohan6.spectrum'];
            this.ompMCoptions.monoEnergy         = 0.1;
            this.ompMCoptions.charge             = 0;
            this.ompMCoptions.sourceGeometry     = 'gaussian';
            this.ompMCoptions.sourceGaussianWidth = 0.1*sigmaGauss;

            % start MC transport
            this.ompMCoptions.dataFolder   = [this.omcFolder filesep 'data' filesep];
            this.ompMCoptions.pegsFile     = [this.omcFolder filesep 'pegs4' filesep '700icru.pegs4dat'];
            this.ompMCoptions.pgs4formFile = [this.omcFolder filesep 'pegs4' filesep 'pgs4form.dat'];

            this.ompMCoptions.global_ecut = 0.7;
            this.ompMCoptions.global_pcut = 0.010;

            % Relative Threshold for dose
            this.ompMCoptions.relDoseThreshold = 1 - this.relativeDosimetricCutOff;

            % Output folders
            this.ompMCoptions.outputFolder = [this.omcFolder filesep 'output' filesep];

            % Create Material Density Cube
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

            % conversion from HU to densities & materials
            for s = 1:dij.numOfScenarios

                HUcube{s} =  matRad_interp3(dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,ct.cubeHU{s}, ...
                    dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest');

                % projecting out of bounds HU values where necessary
                if max(HUcube{s}(:)) > material{end,3}
                    matRad_cfg.dispWarning('Projecting out of range HU values');
                    HUcube{s}(HUcube{s}(:) > material{end,3}) = material{end,3};
                end
                if min(HUcube{s}(:)) < material{1,2}
                    matRad_cfg.dispWarning('Projecting out of range HU values');
                    HUcube{s}(HUcube{s}(:) < material{1,2}) = material{1,2};
                end

                % find material index
                this.cubeMatIx{s} = NaN*ones(dij.doseGrid.dimensions,'int32');
                for i = size(material,1):-1:1
                    this.cubeMatIx{s}(HUcube{s} <= material{i,3}) = i;
                end

                % create an artificial HU lookup table
                hlut = [];
                for i = 1:size(material,1)
                    hlut = [hlut;material{i,2} material{i,4};material{i,3}-1e-10 material{i,5}]; % add eps for interpolation
                end

                this.cubeRho{s} = interp1(hlut(:,1),hlut(:,2),HUcube{s});

            end

            this.ompMCgeo.material = material;

            scale = 10; % to convert to cm

            this.ompMCgeo.xBounds = (dij.doseGrid.resolution.y * (0.5 + [0:dij.doseGrid.dimensions(1)])) ./ scale;
            this.ompMCgeo.yBounds = (dij.doseGrid.resolution.x * (0.5 + [0:dij.doseGrid.dimensions(2)])) ./ scale;
            this.ompMCgeo.zBounds = (dij.doseGrid.resolution.z * (0.5 + [0:dij.doseGrid.dimensions(3)])) ./ scale;

            %% debug visualization
            if this.visBool

                figure
                hold on

                axis equal

                % ct box
                ctCorner1 = [this.ompMCgeo.xBounds(1) this.ompMCgeo.yBounds(1) this.ompMCgeo.zBounds(1)];
                ctCorner2 = [this.ompMCgeo.xBounds(end) this.ompMCgeo.yBounds(end) this.ompMCgeo.zBounds(end)];
                plot3([ctCorner1(1) ctCorner2(1)],[ctCorner1(2) ctCorner1(2)],[ctCorner1(3) ctCorner1(3)],'k' )
                plot3([ctCorner1(1) ctCorner2(1)],[ctCorner2(2) ctCorner2(2)],[ctCorner1(3) ctCorner1(3)],'k' )
                plot3([ctCorner1(1) ctCorner1(1)],[ctCorner1(2) ctCorner2(2)],[ctCorner1(3) ctCorner1(3)],'k' )
                plot3([ctCorner2(1) ctCorner2(1)],[ctCorner1(2) ctCorner2(2)],[ctCorner1(3) ctCorner1(3)],'k' )
                plot3([ctCorner1(1) ctCorner2(1)],[ctCorner1(2) ctCorner1(2)],[ctCorner2(3) ctCorner2(3)],'k' )
                plot3([ctCorner1(1) ctCorner2(1)],[ctCorner2(2) ctCorner2(2)],[ctCorner2(3) ctCorner2(3)],'k' )
                plot3([ctCorner1(1) ctCorner1(1)],[ctCorner1(2) ctCorner2(2)],[ctCorner2(3) ctCorner2(3)],'k' )
                plot3([ctCorner2(1) ctCorner2(1)],[ctCorner1(2) ctCorner2(2)],[ctCorner2(3) ctCorner2(3)],'k' )
                plot3([ctCorner1(1) ctCorner1(1)],[ctCorner1(2) ctCorner1(2)],[ctCorner1(3) ctCorner2(3)],'k' )
                plot3([ctCorner2(1) ctCorner2(1)],[ctCorner1(2) ctCorner1(2)],[ctCorner1(3) ctCorner2(3)],'k' )
                plot3([ctCorner1(1) ctCorner1(1)],[ctCorner2(2) ctCorner2(2)],[ctCorner1(3) ctCorner2(3)],'k' )
                plot3([ctCorner2(1) ctCorner2(1)],[ctCorner2(2) ctCorner2(2)],[ctCorner1(3) ctCorner2(3)],'k' )

                xlabel('x [cm]')
                ylabel('y [cm]')
                zlabel('z [cm]')

                rotate3d on

            end

            %% Create beamlet source
            this.useCornersSCD = true; %false -> use ISO corners

            numOfBixels = [stf(:).numOfRays];
            beamSource = zeros(dij.numOfBeams, 3);

            bixelCorner = zeros(dij.totalNumOfBixels,3);
            bixelSide1 = zeros(dij.totalNumOfBixels,3);
            bixelSide2 = zeros(dij.totalNumOfBixels,3);

            counter = 0;

            for i = 1:dij.numOfBeams % loop over all beams

                % define beam source in physical coordinate system in cm
                beamSource(i,:) = (stf(i).sourcePoint + stf(i).isoCenter)/10;

                for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!

                    counter = counter + 1;

                    dij.beamNum(counter)  = i;
                    dij.rayNum(counter)   = j;
                    dij.bixelNum(counter) = j;

                    if this.useCornersSCD
                        beamletCorners = stf(i).ray(j).rayCorners_SCD;
                    else
                        beamletCorners = stf(i).ray(j).beamletCornersAtIso;
                    end

                    % get bixel corner and delimiting vectors.
                    % a) change coordinate system (Isocenter cs-> physical cs) and units mm -> cm
                    currCorner = (beamletCorners(1,:) + stf(i).isoCenter) ./ scale;
                    bixelCorner(counter,:) = currCorner;
                    bixelSide1(counter,:) = (beamletCorners(2,:) + stf(i).isoCenter) ./ scale - currCorner;
                    bixelSide2(counter,:) = (beamletCorners(4,:) + stf(i).isoCenter) ./ scale - currCorner;

                    if this.visBool
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

            this.ompMCsource.nBeams = dij.numOfBeams;
            this.ompMCsource.iBeam = dij.beamNum(:);

            % Switch x and y directions to match ompMC cs.
            this.ompMCsource.xSource = beamSource(:,2);
            this.ompMCsource.ySource = beamSource(:,1);
            this.ompMCsource.zSource = beamSource(:,3);

            this.ompMCsource.nBixels = sum(numOfBixels(:));
            this.ompMCsource.xCorner = bixelCorner(:,2);
            this.ompMCsource.yCorner = bixelCorner(:,1);
            this.ompMCsource.zCorner = bixelCorner(:,3);

            this.ompMCsource.xSide1 = bixelSide1(:,2);
            this.ompMCsource.ySide1 = bixelSide1(:,1);
            this.ompMCsource.zSide1 = bixelSide1(:,3);

            this.ompMCsource.xSide2 = bixelSide2(:,2);
            this.ompMCsource.ySide2 = bixelSide2(:,1);
            this.ompMCsource.zSide2 = bixelSide2(:,3);

            if this.visBool
                plot3(this.ompMCsource.ySource,this.ompMCsource.xSource,this.ompMCsource.zSource,'rx')
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



    end


end

