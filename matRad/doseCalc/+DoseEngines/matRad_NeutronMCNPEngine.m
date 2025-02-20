classdef matRad_NeutronMCNPEngine < DoseEngines.matRad_MonteCarloEngineAbstract
    % Engine for neutron dose calculation using monte carlo calculation
    % specificly MCNP
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
        possibleRadiationModes = {'neutrons'};
        name = 'MCNP dose engine';
        shortName = 'MCNP';
    end

    properties
        config;             %Holds an instance of all configurable parameters

        %Other Dose Calculation Properties
        externalCalculation = true;
        
        % Question dialogs avoided
        useDICOMinfoRescale = true;
        useLungQuestionDialog = false;
        
        calcRMFparameters = false;
        MCNPinstallationCheck;
        MCNPFolder;
    end

    properties (SetAccess = protected, GetAccess = public)

        currFolder = pwd; %folder path when set

        constantRBE = NaN;              % constant RBE value
    end

    methods

        function this = matRad_NeutronMCNPEngine(pln)
            % Constructor
            %
            % call
            %   engine = DoseEngines.matRad_NeutronMCNPEngine(ct,stf,pln,cst)
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

            this.config = matRad_MCNPConfig();

            % check if bio optimization is needed and set the
            % coresponding boolean accordingly
            % TODO:
            % This should not be handled here as an optimization property
            % We should rather make optimization dependent on what we have
            % decided to calculate here.
            if nargin > 0
                if (isfield(pln,'propOpt')&& isfield(pln.propOpt,'bioOptimization')&& ...
                        isequal(pln.propOpt.bioOptimization,'RBExSecPartDose_MCDS_RMFmodel'))
                    this.calcBioDose = true;
                elseif strcmp(pln.radiationMode,'neutrons') && isfield(pln,'propOpt') && isfield(pln.propOpt,'bioOptimization') && isequal(pln.propOpt.bioOptimization,'const_RBExD')
                    this.constantRBE = 4;
                end
            end
        end

        function setDefaults(this)
            this.setDefaults@DoseEngines.matRad_MonteCarloEngineAbstract();

            % future code for property validation on creation here
            matRad_cfg = MatRad_Config.instance();

            %Assign default parameters from MatRad_Config
            %this.doseGrid.resolution    = matRad_cfg.propDoseCalc.defaultResolution;
            %this.multScen = 'nomScen';
            
            %Set Default MCNP path
            %Set folder
            this.MCNPFolder = [matRad_cfg.matRadRoot filesep 'MCNP'];        
        end
    end

    methods(Access = protected)

        function dij = calcDose(this,ct,cst,stf)
            % matRad MCNP monte carlo dose calculation wrapper
            %
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
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            matRad_cfg = MatRad_Config.instance();

            %Now we can run initDoseCalc as usual
            dij = this.initDoseCalc(ct,cst,stf);

            % For MCNP TMESH calculations the ct grid has to be downsampled to the dose grid
            for s = 1:dij.numOfScenarios
                if s>1
                    error('Multiple scenarios not yet implemented for MCNP dose calculations.')
                end
                ct.doseGridCT.HUcube{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cubeHU{s}, ...
                    dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
                cst_doseGrid = matRad_resizeCstToGrid(cst,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,...
                    dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);
            end

            % set absolute calibration factor
            % Convert MeV/g to J/kg, output is now in Gy/source particle
            absCalibrationFactorMCNP = 1.602177e-19*1e6*1e3;

            %matRad_cfg.matRadRoot 

            % Generate log-file
            pathLog = strcat(this.MCNPFolder, filesep, 'logFile');

            try diary(fullfile(pathLog, strcat(matRad_getTime4log, '_neutronDoseCalculation')))
            catch
                mkdir('logFile')
                diary(fullfile(pathLog, strcat(matRad_getTime4log, '_neutronDoseCalculation')))
            end

            % Load predefined conversion properties for tissue characterization
            % according to CT values - elemental composition will be assigned according
            % to these predefined tissue intervals later
            matRad_cfg.dispInfo('*****\n')
            matRad_cfg.dispInfo('Load pre-defined HU conversion properties and MCNP cross sections from conversionCT2tissue.mat.\n')
            matRad_cfg.dispInfo('Note: Modification of conversionCT2tissue.mat using matRad_genVarCT2tissueConversion.m.\n')
            load([this.MCNPFolder,filesep,'conversionCT2tissue.mat'])
            matRad_cfg.dispInfo('*****\n')

            %% Process CT data
            % Check ct for MCNP Simulation
            if ct.resolution.x ~= ct.resolution.y
                error('x- and y-resolution have to be equal for the simulation.')
            end

            % Set HU outside body to air, i.e. neglect everything outside body for the
            % simulation
            % Note: body structure is the only normal tissue structure that
            % has to be contured

            bodyStructureName = 'Body'; % Default name for body structure is 'Body'
            try
                cstBodyIndex = matRad_findBodyStructureCST(cst, bodyStructureName);
            catch
                prompt = {'Please enter body structure name:'};
                dlgtitle = 'Find Body Structure';
                bodyStructureName = inputdlg(prompt,dlgtitle);
                cstBodyIndex = matRad_findBodyStructureCST(cst, bodyStructureName);
            end

            % Process HU values
            matRad_cfg.dispInfo('*****\n')
            matRad_cfg.dispInfo(['Properties from (scaled) HU loaded are: Minimum value: ',num2str(min(ct.cubeHU{1}, [], 'all')), ' and \n'])
            matRad_cfg.dispInfo(['Maximum value: ',num2str(max(ct.cubeHU{1}(cst{cstBodyIndex,4}{1}), [], 'all')), '.\n'])
            if ~isfield(ct, 'dicomInfo') || ~isfield(ct.dicomInfo, 'RescaleSlope') || ~isfield(ct.dicomInfo, 'RescaleIntercept')
                matRad_cfg.dispWarning('No information on rescale slope and/or intercept provided in DICOM data. Calculation might crash...\n')
            else
                matRad_cfg.dispInfo(['Rescale slope from CT data is read to be: ',num2str(ct.dicomInfo.RescaleSlope), ' and rescale intercept is read to be ',num2str(ct.dicomInfo.RescaleIntercept),'.\n'])
            end

            % Set values outside body to air
            try
                cstTargetIndex = matRad_findTargetStructureCST(cst);
            catch
                matRad_cfg.dispError('Target structure has to be set in matRad.\n')
            end

            maskNonBody = ones(ct.cubeDim);
            maskNonBody_doseGrid = ones(dij.doseGrid.dimensions);
            bodyIdx = [cst{sort([cstTargetIndex, cstBodyIndex]),4}];
            bodyIdx = unique(vertcat(bodyIdx{:}));
            bodyIdx_doseGrid = [cst_doseGrid{sort([cstTargetIndex, cstBodyIndex]),4}];
            bodyIdx_doseGrid = unique(vertcat(bodyIdx_doseGrid{:}));

            maskNonBody(bodyIdx) = 0;
            maskNonBody_doseGrid(bodyIdx_doseGrid) = 0;
            ct.cube{1}(maskNonBody>0) = 0;
            ct.cubeHU{1}(maskNonBody>0) = 0;
            ct.doseGridCT.HUcube{1}(maskNonBody_doseGrid>0) = 0;


            if ~isprop(this, 'useDICOMinfoRescale')
                matRad_cfg.dispInfo('Please use question dialog to decide how to convert to scaled HU.\n')
                matRad_cfg.dispInfo('*****\n')
                % DICOM rescaling
                answer = questdlg('Would you like to use DICOM rescale slope and intercept? If not, an offset of 1000 will be added to the HU values to get re-scale HUs.', ...
                    'Use DICOM Info', ...
                    'Yes', 'No', 'No');

                switch answer
                    case 'Yes'
                        matRad_cfg.dispInfo('*****\n')
                        matRad_cfg.dispInfo('You decided to use the following parameters to re-scale (scaled) HU data given in ct.cubeHU to HU in ct.cube.\n')
                        matRad_cfg.dispInfo(['Rescale HU: slope=', num2str(ct.dicomInfo.RescaleSlope), ' intercept=', num2str(ct.dicomInfo.RescaleIntercept),'\n'])
                        matRad_cfg.dispInfo('*****\n')
                        ct.doseGridCT.HUcube{1} = ct.doseGridCT.HUcube{1}.*ct.dicomInfo.RescaleSlope + abs(ct.dicomInfo.RescaleIntercept);
                    case 'No'
                        matRad_cfg.dispInfo('*****\n')
                        matRad_cfg.dispInfo('You decided not to use DICOM rescale slope and intercept.\n')
                        matRad_cfg.dispInfo('*****\n')
                        ct.doseGridCT.HUcube{1} = ct.doseGridCT.HUcube{1} + 1000;
                end
                clear answer
            elseif this.useDICOMinfoRescale
                matRad_cfg.dispInfo('*****\n')
                matRad_cfg.dispInfo('Re-scaling of (scaled) HU data given in ct.cubeHU to HU in ct.cube selected according to property useDICOMinfoRescale.\n')
                matRad_cfg.dispInfo('DICOM rescale slope and intercept are used to re-scale.\n')
                matRad_cfg.dispInfo(['Rescale HU: slope=', num2str(ct.dicomInfo.RescaleSlope), ' intercept=', num2str(ct.dicomInfo.RescaleIntercept),'\n'])
                matRad_cfg.dispInfo('*****\n')
                ct.doseGridCT.HUcube{1} = ct.doseGridCT.HUcube{1}.*ct.dicomInfo.RescaleSlope + abs(ct.dicomInfo.RescaleIntercept);
            elseif ~this.useDICOMinfoRescale
                matRad_cfg.dispInfo('*****\n')
                matRad_cfg.dispInfo('Re-scaling of (scaled) HU data given in ct.cubeHU to HU in ct.cube selected according to property useDICOMinfoRescale.\n')
                matRad_cfg.dispInfo('You decided not to use DICOM rescale slope and intercept. An offset of +1000 is added.\n')
                matRad_cfg.dispInfo('*****\n')
                ct.doseGridCT.HUcube{1} = ct.doseGridCT.HUcube{1} + 1000;

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Additional Information on Conversion                                                                        %
            % Re-scale HU from intensity values given in DICOM by using 'rescale intercept' and 'rescale slope'           %
            % Def. HU: HU = 1000*((mu - mu_water)/mu_water) -> HU_water = 0                                               %
            % Note: For our purpose, HU values should start at zero.                                                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Infere tissue characteristics from HU
            matRad_cfg.dispInfo('*****\n')
            matRad_cfg.dispInfo(['Properties from rescaled HU are: Minimum value: ',num2str(min(ct.doseGridCT.HUcube{1}, [], 'all')), ' and \n'])
            matRad_cfg.dispInfo(['Maximum value: ',num2str(max(ct.doseGridCT.HUcube{1}, [], 'all')), '.\n'])
            matRad_cfg.dispInfo('*****\n')
            matRad_cfg.dispInfo('Material types are assigned using the following HU intervals:\n')
            maxHUbin_nonEmpty = 1;  % Find last non-empty entry
            while ~isempty(binIntervals(maxHUbin_nonEmpty+1).HUbin) && (maxHUbin_nonEmpty <= size(binIntervals,2))
                maxHUbin_nonEmpty = maxHUbin_nonEmpty +1;
            end
            for i=1:maxHUbin_nonEmpty
                matRad_cfg.dispInfo([binIntervals(i).name, ': ', num2str(binIntervals(i).HUbin(1)), ' to ', num2str(binIntervals(i).HUbin(2)),'\n']);
            end
            matRad_cfg.dispInfo('*****\n')

            [cst_doseGrid, ct.doseGridCT.tissueBin] = matRad_segmentationCTscan(ct.doseGridCT.HUcube{1}, dij.doseGrid.resolution, binIntervals, cst_doseGrid, cstBodyIndex, cstTargetIndex, this.useLungQuestionDialog);
            ct.doseGridCT.cubeDim = dij.doseGrid.dimensions;

            % Calculate density for CT voxels and resize afterwards -> caution: step-
            % wise definition of conversion causes a difference concerning the ordering
            % of the conversion and resize operation
            matRad_cfg.dispInfo('*****\n')
            matRad_cfg.dispInfo('Calculate density from CT data with density given in [g/cm^3]\n')
            matRad_cfg.dispInfo('*****\n')
            ct.doseGridCT.density{1} = hounsfield2density(ct.doseGridCT.HUcube{1})*1e-3; % rescale from kg/m^3 to g/cm^3
            ct.doseGridCT.density{1}(ct.doseGridCT.density{1}<=(hounsfield2density(segVar.upperLimitAir)*1e-3)) = segVar.densityAir;

            %% Rescale ct information to cm for MCNP runfile
            varHelper.rescaleFactor = 1e-1; % conversion from mm to cm
            ct.doseGridCT.x_MCNP = dij.doseGrid.resolution.x*varHelper.rescaleFactor;
            ct.doseGridCT.y_MCNP = dij.doseGrid.resolution.y*varHelper.rescaleFactor;
            ct.doseGridCT.z_MCNP = dij.doseGrid.resolution.z*varHelper.rescaleFactor;


            %% Create MCNP runfile blocks A and B
            varHelper.simPropMCNP.loopCounter = false; % try to generate MCNP runfile with maximum 99999 elements for reasons of performance
            varHelper.simPropMCNP.MCNP_limitNumberOfElements = 99999-1; % minus one since we need one integer for the source surface

            [varHelper.simPropMCNP.control_makeTargetMCNP, varHelper.simPropMCNP.fileID_A, varHelper.simPropMCNP.fileID_B, varHelper.simPropMCNP.geometryOption] = ...
                matRad_makeTargetMCNP(ct, varHelper.simPropMCNP);

            if ~varHelper.simPropMCNP.control_makeTargetMCNP
                varHelper.simPropMCNP.MCNP_limitNumberOfElements = 99999999-1;  % minus one since we need one integer for the source surface
                varHelper.simPropMCNP.loopCounter = true;
                [varHelper.simPropMCNP.control_makeTargetMCNP,~ ,~ ,~ ] = ...
                    matRad_makeTargetMCNP(ct, varHelper.simPropMCNP);
            end

            if ~varHelper.simPropMCNP.control_makeTargetMCNP
                matRad_cfg.dispError('Number of defined elements for simulation is too high! MCNP6 only allows 99,999,999 cells and surfaces in total.\n')
            end

            %% Create MCNP runfile block C (source etc.)
            % Set default number of particles
            varHelper.simPropMCNP.numberParticles = this.config.Num_Primaries;

            % Get total number of bixels and write source card
            % Total number of bixel is counter for j in dij matrix

            % C.1 Source
            varHelper.totalNumberBixels = 0;
            for counterField =1:size(stf,2)
                varHelper.simPropMCNP.counterField = counterField;
                for counterRay=1:stf(counterField).numOfRays
                    varHelper.simPropMCNP.counterRay = counterRay;
                    % Calculate position of MLC field in LPS system
                    stf(counterField).ray(counterRay).rayPosMLC = stf(counterField).isoCenter + stf(counterField).ray(counterRay).rayPos + stf(counterField).sourcePoint;
                    varHelper.totalNumberBixels = varHelper.totalNumberBixels + 1;

                    % Generate source card for each bixel - see also description of
                    % matRad^_makeSourceMCNP(...)
                    [control_makeSourceMCNP, varHelper] = matRad_makeSourceMCNP(this, stf, varHelper, counterField, counterRay);
                end
            end

            varHelper.simPropMCNP.sourceBlockNames = strings(1,varHelper.totalNumberBixels);
            for cntList=1:varHelper.totalNumberBixels; varHelper.simPropMCNP.sourceBlockNames(cntList)=['blockC_source', int2str(cntList)]; end

            % Open file to write rest of block C
            % Source and the rest of block C input are seperated so the source
            % positioning can be done easily w/o wasting time on redundant writing of
            % the rest (like MODE and PHYS card) into a text file.

            pathRunfiles = strcat(matRad_cfg.matRadRoot,filesep, 'MCNP', filesep, 'runfiles_tmp', filesep);
            fileID_C_rest = fopen(strcat(pathRunfiles,'blockC_rest'), 'w');

            % C.2 Physics and problem termination
            fprintf(fileID_C_rest, 'C ***************************************************************\n');
            fprintf(fileID_C_rest, 'C C.2: Physics\n');
            fprintf(fileID_C_rest, 'C ***************************************************************\n');

            matRad_definePhysicsMCNP(fileID_C_rest, this, binIntervals, varHelper.simPropMCNP); % Define PHYS-card

            % C.3 Materials
            fprintf(fileID_C_rest, 'C ***************************************************************\n');
            fprintf(fileID_C_rest, 'C C.3: Materials\n');
            fprintf(fileID_C_rest, 'C ***************************************************************\n');

            matIdentifierTxt = 'M%d $ %s \n';
            matCompositionTxt = '           %d%s %.9f \n';

            for counterMaterial = 1:size(binIntervals,2)
                fprintf(fileID_C_rest, matIdentifierTxt, counterMaterial, binIntervals(counterMaterial).name);
                fprintf(fileID_C_rest,'	   plib=14p\n');
                fprintf(fileID_C_rest,'	   hlib=70h\n');
                for counterComponent = 1:size(binIntervals(counterMaterial).ZAID,2)
                    fprintf(fileID_C_rest,matCompositionTxt, binIntervals(counterMaterial).ZAID(counterComponent), crossSectionsLibrary(binIntervals(counterMaterial).crossSection(counterComponent),:), binIntervals(counterMaterial).percentageMass(counterComponent));
                end
            end
            clear counterMaterial; clear counterComponent;

            % C.4 Tally
            matRad_makeTallyMCNP(this, ct, fileID_C_rest, binIntervals)

            fclose(fileID_C_rest);

            %% Concatenate all blocks to one runfile for each ray
            matRad_concatenateRunfiles(varHelper, pathRunfiles);

            %% Run MCNP calculation
            matRad_bixelDoseCalculatorMCNP(this);

            %% Evaluate MCNP results
            dij = matRad_evaluateTallyMCNP(dij, cst, ct);
            
            %% Switch off diary
            diary off
        end

        function setBinaries(this)
            % setBinaries check if the binaries are available on the current
            % machine and sets to the mcsquarebinary object property
            %

            binaryFound = this.checkBinaries();
            this.MCNPinstallationCheck = binaryFound;
        end

        function dij = initDoseCalc(this,ct,cst,stf)
            %% Assingn and check parameters
            matRad_cfg = MatRad_Config.instance();

            % check if binaries are available
            % Executables for simulation
            this.setBinaries();

            % set and change to MCsquare binary folder
            this.currFolder = pwd;

            %% Call Superclass init function
            dij = initDoseCalc@DoseEngines.matRad_MonteCarloEngineAbstract(this,ct,cst,stf);

            % Since MCNP setup only allows similar resolution in x&y, we do some
            % extra checks on that before calling the normal initDoseCalc. 
            if dij.doseGrid.resolution.x ~= dij.doseGrid.resolution.y
                matRad_cfg.dispError('Voxel size in x and y do not agree.\n');
            end

            %% Validate and preset some additional dij variables
            %Issue a warning when we have more than 1 scenario
            if dij.numOfScenarios ~= 1
                matRad_cfg.dispWarning('MCNP is only implemented for single scenario use at the moment. Will only use the first Scenario for Monte Carlo calculation!');
            end

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
        function binaryFound = checkBinaries()
            % checkBinaries check if MCNP is installed on the machine and
            % path variables are properly set to run MCNP in matRad
            matRad_cfg = MatRad_Config.instance();                       

            if ispc
                [~,cmdout] = system('mcnp6');
                if strcmp(cmdout(2:5), 'mcnp')
                    binaryFound = true;
                else
                    matRad_cfg.dispWarning('Could not test MCNP. Please check installation and path variables.\n');
                end
            elseif ismac
                matRad_cfg.dispWarning('Check for MCNP installation not yet implemented. Check set to false.\n');
                binaryFound = false;
            elseif isunix
                matRad_cfg.dispWarning('Check for MCNP installation not yet implemented. Check set to false.\n');
                binaryFound = false;

            else
                binaryFound = false;
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
                checkBasic = any(strcmp(machine.meta.calculationMode, 'MonteCarlo'));

                %check modality
                checkModality = any(strcmp(DoseEngines.matRad_NeutronMCNPEngine.possibleRadiationModes, machine.meta.radiationMode));

                preCheck = checkBasic && checkModality;

                if ~preCheck
                    return;
                end
                available = true;
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return;
            end

        end

    end
end

