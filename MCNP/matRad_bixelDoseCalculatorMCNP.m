function dij = matRad_bixelDoseCalculatorMCNP(this)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad neutron dose calculation for each bixel individually
%
% Neutron dose engine A: Monte Carlo - MCNP6
%
% call
%   dij = matRad_calcPhotonDose(pathToRunfiles, stf, ct, pln, cst, binIntervals)
%
% input
%   pathToRunfiles: indicate path to MCNP runfiles here
%   stf, ct, pln, binIntervals
%
% output
%   dij:            matRad dij struct
%
% References
%   [1] PELOWITZ, D. B., et al. MCNP6 Users Manual. LACP-00634, May, 2013.
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 11/2018
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matRad_cfg = MatRad_Config.instance();

if this.MCNPinstallationCheck && ~this.externalCalculation
    %% Go to runfiles and get list of runfiles within directory
    oldPath = pwd;
    cd(pathToRunfiles)
    runFileList = dir('MCNPrunfile_bixel*');

    wb = waitbar(0, ['Calculating dose for bixel: ', num2str(1)], 'Name', 'Dose Calculation with MCNP');
    numberCores4U = feature('numcores');    %get number of physical cores available on machine

    %% Run calculation for every bixel
    for bixelCounter=1:size(runFileList,1)

        disp('*****')
        disp(['MCNP calculation of dose distribution for bixel ', num2str(bixelCounter), '...'])
        disp('*****')

        tic;

        %     waitbar(bixelCounter/size(runFileList,1), wb, ['Calculating dose for bixel: ', num2str(bixelCounter)], 'Name', 'Dose Calculation with MCNP');
        if ispc % Set up MCNP environment and run MCNP in one command, make sure mcnp_env_620.bat from the installation exist in the runfile directiory
            system(['mcnp_env_620.bat & mcnp6 I=', runFileList(bixelCounter).name, ...
                ' OUTP=', runFileList(bixelCounter).name, 'o ', ...
                ' RUNTPE=', runFileList(bixelCounter).name, 'r ', ...
                ' MCTAL=', runFileList(bixelCounter).name, 'm ', ...
                ' tasks ' ,num2str(numberCores4U)])
            %             ' SRCTP=', runFileList(bixelCounter).name, 's ', ...
            %             ' PTRAC=', runFileList(bixelCounter).name, 'p ', ...
            %             ' MESHTAL=', runFileList(bixelCounter).name, 'meshtal ', ...

            % Clean up
            delete(strcat(runFileList(bixelCounter).name, 'o'))
            delete(strcat(runFileList(bixelCounter).name, 'r'));
            %         delete(strcat(runFileList(bixelCounter).name, 's'))
            %         delete(strcat(runFileList(bixelCounter).name, 'p'))

        else    % Usually MCNP environment is set up permanently e.g. on linux systems, if not do so or change next system command
            system(['mcnp6 I=', runFileList(bixelCounter).name, ...
                ' OUTP=', runFileList(bixelCounter).name, 'o ', ...
                ' RUNTPE=', runFileList(bixelCounter).name, 'r ', ...
                ' SRCTP=', runFileList(bixelCounter).name, 's ', ...
                ' MCTAL=', runFileList(bixelCounter).name, 'm ', ...
                ' PTRAC=', runFileList(bixelCounter).name, 'p ', ...
                ' MESHTAL=', runFileList(bixelCounter).name, 'meshtal ', ...
                ' tasks ' ,num2str(numberCores4U)])
            % Clean up
            delete(strcat(runFileList(bixelCounter).name, 'o'))
            delete(strcat(runFileList(bixelCounter).name, 'r'));
            delete(strcat(runFileList(bixelCounter).name, 's'))
            delete(strcat(runFileList(bixelCounter).name, 'p'))
        end

        calculationTime = toc;
        disp('*****')
        disp(['Calculation for bixel ', num2str(bixelCounter), ' took ', num2str(calculationTime), ' seconds.'])
        disp('*****')

    end

    close(wb)
elseif this.externalCalculation
    matRad_cfg.dispInfo('Please use question dialog to continue after finishing external calculations.\n')
    matRad_cfg.dispInfo('*****\n')
    % External calculation
    answer = questdlg('Did external MCNP simulations finish?', ...
        'External Calculation', ...
        'Yes', 'No', 'No');

    while ~strcmp(answer, 'Yes')
        matRad_cfg.dispInfo('matRad will crash if you continue without finishing external calculation.\n')
        answer = questdlg('Did external MCNP simulations finish?', ...
            'External Calculation', ...
            'Yes', 'No', 'No');
    end
elseif ~this.externalCalculation && this.MCNPinstallationCheck
    matRad_cfg.dispWarning('MCNP simulation requested but no MCNP installation found on your computer!\n')
end
%% Read output from MCNP calculation and generate dij matix
if strcmp(pln.propMCNP.tallySpecifier, 'KERMA_F4')
    % Get list of meshtally data
    tallyDataList = dir('MCNPrunfile_bixel*meshtal');
    % Re-organize list
    dummyList = struct;
    lengthList = sum([stf(:).numOfRays]);
    for listCounter = 1:lengthList
        i=1; while i<= lengthList && ~strcmp(tallyDataList(i).name, strcat('MCNPrunfile_bixel', int2str(listCounter), 'meshtal')); i=i+1; end
        dummyList(listCounter).name = tallyDataList(i).name;
        dummyList(listCounter).folder = tallyDataList(i).folder;
        dummyList(listCounter).date = tallyDataList(i).date;
        dummyList(listCounter).bytes = tallyDataList(i).bytes;
        dummyList(listCounter).isdir = tallyDataList(i).isdir;
        dummyList(listCounter).datenum = tallyDataList(i).datenum;
    end
    tallyDataList = dummyList; clear dummyList;
elseif strcmp(pln.propMCNP.tallySpecifier, 'TotalDose_TMESH')
    % Get list of mctal data
    tallyDataList = dir('MCNPrunfile_bixel*m');
    % Re-organize list
    dummyList = struct;
    lengthList = sum([stf(:).numOfRays]);
    for listCounter = 1:lengthList
        i=1; while i<= lengthList && ~strcmp(tallyDataList(i).name, strcat('MCNPrunfile_bixel', int2str(listCounter), 'm')); i=i+1; end
        dummyList(listCounter).name = tallyDataList(i).name;
        dummyList(listCounter).folder = tallyDataList(i).folder;
        dummyList(listCounter).date = tallyDataList(i).date;
        dummyList(listCounter).bytes = tallyDataList(i).bytes;
        dummyList(listCounter).isdir = tallyDataList(i).isdir;
        dummyList(listCounter).datenum = tallyDataList(i).datenum;
    end
end

% Set meta information for dij
dij.numOfBeams         = pln.propStf.numOfBeams;
dij.numOfVoxels        = prod(ct.cubeDim);
dij.resolution         = ct.resolution;
dij.dimensions         = ct.cubeDim;
dij.numOfScenarios     = 1;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);

dij.bixelNum = NaN*ones(dij.totalNumOfBixels,1); % dij.totalNumOfBixels
dij.rayNum   = NaN*ones(dij.totalNumOfBixels,1); % equals the number of
dij.beamNum  = NaN*ones(dij.totalNumOfBixels,1); % columnes for the dij
% matrix

% Allocate space for dij.physicalDose sparse matrix
% for i = 1:dij.numOfScenarios   % more than one scenario not implemented
% in matRad
%     dij.physicalDose{i} = spalloc(prod(ct.cubeDim),dij.totalNumOfBixels,1);
% end

% Allocate dij dose matrix
dij.physicalDose{1} = spalloc(prod(ct.cubeDim),dij.totalNumOfBixels,1);
for counterBixel=1:dij.totalNumOfBixels
    dij.doseMatrixBixel(counterBixel).photonDose{1} = spalloc(prod(ct.cubeDim),1,1);
    dij.doseMatrixBixel(counterBixel).neutronDose{1} = spalloc(prod(ct.cubeDim),1,1);
end

if strcmp(pln.propMCNP.tallySpecifier, 'KERMA_F4')
    counterDijColumns = 0; % counter to go through the columns of dij

    for counterBeam = 1:dij.numOfBeams
        for counterBixel = 1:dij.numOfRaysPerBeam(counterBeam)
            counterDijColumns = counterDijColumns+1;

            dij.bixelNum(counterDijColumns) = counterBixel;
            dij.rayNum(counterDijColumns) = counterBixel;
            dij.beamNum(counterDijColumns) = counterBeam;

            % Evaluate mesh tally and store neutron and photon KERMA
            disp('*****')
            disp('Evaluate Tally Data.')
            disp('*****')
            tic;
            doseMatrixBixel = matRad_evaluateMeshTallyMCNP(stf, tallyDataList(counterDijColumns).name, pln, binIntervals, ct);
            disp(['Reading out tally data took: ', num2str(toc), ' seconds.'])
            disp('*****')

            dij.doseMatrixBixel(counterDijColumns).photonDose{1} = ...
                sparse(find(doseMatrixBixel.photonDose),1,doseMatrixBixel.photonDose(doseMatrixBixel.photonDose~=0),dij.numOfVoxels,1);
            dij.doseMatrixBixel(counterDijColumns).photonDose_relError{1} = ...
                sparse(find(doseMatrixBixel.photonDose_relError),1,doseMatrixBixel.photonDose_relError(doseMatrixBixel.photonDose_relError~=0),dij.numOfVoxels,1);
            dij.doseMatrixBixel(counterDijColumns).neutronDose{1} = ...
                sparse(find(doseMatrixBixel.neutronDose),1,doseMatrixBixel.neutronDose(doseMatrixBixel.neutronDose~=0),dij.numOfVoxels,1);
            dij.doseMatrixBixel(counterDijColumns).neutronDose_relError{1} = ...
                sparse(find(doseMatrixBixel.neutronDose_relError),1,doseMatrixBixel.neutronDose_relError(doseMatrixBixel.neutronDose_relError~=0),dij.numOfVoxels,1);

            % Store total KERMA from neutron and photon interaction
            doseContainer = doseMatrixBixel.photonDose + doseMatrixBixel.neutronDose;
            dij.physicalDose{1}(:,counterDijColumns) = sparse(find(doseContainer), 1, doseContainer(doseContainer~=0), dij.numOfVoxels,1);

            % Calculate relative error for sum
            doseContainer_relError = sqrt((doseMatrixBixel.photonDose.*doseMatrixBixel.photonDose_relError).^2 + ...
                (doseMatrixBixel.neutronDose.*doseMatrixBixel.neutronDose_relError).^2) ./ ...
                doseContainer;
            doseContainer_relError(doseContainer==0)=0; % Division by doseContainer with zero values leads to NaN, where doseContainer is zero, the error is zero too
            dij.physicalDose_relError{1}(:,counterDijColumns) = sparse(find(doseContainer_relError), 1, doseContainer_relError(doseContainer_relError~=0), dij.numOfVoxels,1);

            % Save relative error for each radiotherapy structure and bixel
            for counterRTStruct=1:size(cst,1)
                dij.relError_radioTherpyStruct(counterRTStruct,counterDijColumns).name = cst{counterRTStruct,2};
                dij.relError_radioTherpyStruct(counterRTStruct,counterDijColumns).bixelNumber = counterDijColumns;
                dij.relError_radioTherpyStruct(counterRTStruct,counterDijColumns).meanError = mean(doseContainer_relError(cst{counterRTStruct,4}{1}));
                dij.relError_radioTherpyStruct(counterRTStruct,counterDijColumns).maxError = max(doseContainer_relError(cst{counterRTStruct,4}{1}));
                dij.relError_radioTherpyStruct(counterRTStruct,counterDijColumns).medianError = median(doseContainer_relError(cst{counterRTStruct,4}{1}));
            end

        end
    end
elseif strcmp(pln.propMCNP.tallySpecifier, 'TotalDose_TMESH')
    counterDijColumns = 0; % counter to go through the columns of dij

    for counterBeam = 1:dij.numOfBeams
        for counterBixel = 1:dij.numOfRaysPerBeam(counterBeam)
            counterDijColumns = counterDijColumns+1;

            dij.bixelNum(counterDijColumns) = counterBixel;
            dij.rayNum(counterDijColumns) = counterBixel;
            dij.beamNum(counterDijColumns) = counterBeam;

            % Evaluate mesh tally and store neutron and photon dose
            disp('*****')
            disp('Evaluate Tally Data.')
            disp('*****')
            tic;
            doseMatrixBixel = matRad_evaluateMeshTallyMCNP(stf, tallyDataList(counterDijColumns).name, pln, binIntervals, ct);
            disp(['Reading out tally data took: ', num2str(toc), ' seconds.'])
            disp('*****')

            dij.doseMatrixBixel(counterDijColumns).photonDose{1} = ...
                sparse(find(doseMatrixBixel.photonDose),1,doseMatrixBixel.photonDose(doseMatrixBixel.photonDose~=0),dij.numOfVoxels,1);
            dij.doseMatrixBixel(counterDijColumns).photonDose_relError{1} = ...
                sparse(find(doseMatrixBixel.photonDose_relError),1,doseMatrixBixel.photonDose_relError(doseMatrixBixel.photonDose_relError~=0),dij.numOfVoxels,1);
            dij.doseMatrixBixel(counterDijColumns).neutronDose{1} = ...
                sparse(find(doseMatrixBixel.neutronDose),1,doseMatrixBixel.neutronDose(doseMatrixBixel.neutronDose~=0),dij.numOfVoxels,1);
            dij.doseMatrixBixel(counterDijColumns).neutronDose_relError{1} = ...
                sparse(find(doseMatrixBixel.neutronDose_relError),1,doseMatrixBixel.neutronDose_relError(doseMatrixBixel.neutronDose_relError~=0),dij.numOfVoxels,1);
            dij.doseMatrixBixel(counterDijColumns).physicalDose_relError{1} = ...
                sparse(find(doseMatrixBixel.physicalDose_relError),1,doseMatrixBixel.physicalDose_relError(doseMatrixBixel.physicalDose_relError~=0),dij.numOfVoxels,1);

            % Store total dose from TMESH tally
            dij.physicalDose{1}(:,counterDijColumns) = ...
                sparse(find(doseMatrixBixel.physicalDose), 1, doseMatrixBixel.physicalDose(doseMatrixBixel.physicalDose~=0), dij.numOfVoxels,1);

            % Save relative error for each radiotherapy structure and bixel
            for counterRTStruct=1:size(cst,1)
                dij.relError_radioTherpyStruct(counterRTStruct,counterDijColumns).name = cst{counterRTStruct,2};
                dij.relError_radioTherpyStruct(counterRTStruct,counterDijColumns).bixelNumber = counterDijColumns;
                dij.relError_radioTherpyStruct(counterRTStruct,counterDijColumns).meanError = mean(doseMatrixBixel.physicalDose_relError(cst{counterRTStruct,4}{1}));
                dij.relError_radioTherpyStruct(counterRTStruct,counterDijColumns).maxError = max(doseMatrixBixel.physicalDose_relError(cst{counterRTStruct,4}{1}));
                dij.relError_radioTherpyStruct(counterRTStruct,counterDijColumns).medianError = median(doseMatrixBixel.physicalDose_relError(cst{counterRTStruct,4}{1}));
            end

            %% Calculate RMF parameters
            if isfield(pln.propOpt,'bioOptimization')
                if strcmp(pln.propOpt.bioOptimization,'RBExSecPartDose_MCDS_RMFmodel')
                    [dij, cst] = matRad_getRMFmodelParameters4neutrons(dij, counterDijColumns, doseMatrixBixel, pln, ct, cst);
                else
                    %                     error('Biological optimization mode unknown for neutrons.')
                end
            end


        end
    end
end

%%
cd(oldPath)