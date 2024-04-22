function dij = matRad_evaluateTallyMCNP(dij, cst)
%% Read output from MCNP calculation and generate dij matix
% Get list of mctal data
matRad_cfg = MatRad_Config.instance();
cd(strcat(matRad_cfg.matRadRoot, filesep, 'MCNP', filesep, 'runfiles_tmp'));
tallyDataList = dir('MCNPrunfile_bixel*m');
% Re-organize list
dummyList = struct;
lengthList = dij.totalNumOfRays;
for listCounter = 1:lengthList
    i=1; while i<= lengthList && ~strcmp(tallyDataList(i).name, strcat('MCNPrunfile_bixel', int2str(listCounter), 'm')); i=i+1; end
    dummyList(listCounter).name = tallyDataList(i).name;
    dummyList(listCounter).folder = tallyDataList(i).folder;
    dummyList(listCounter).date = tallyDataList(i).date;
    dummyList(listCounter).bytes = tallyDataList(i).bytes;
    dummyList(listCounter).isdir = tallyDataList(i).isdir;
    dummyList(listCounter).datenum = tallyDataList(i).datenum;
end

% Allocate dij dose matrix
dij.physicalDose{1} = spalloc(prod(dij.doseGrid.dimensions),dij.totalNumOfBixels,1);
for counterBixel=1:dij.totalNumOfBixels
    dij.doseMatrixBixel(counterBixel).physicalDose_relError{1} = spalloc(prod(dij.doseGrid.dimensions),1,1);
end

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
        doseMatrixBixel = matRad_evaluateMeshTallyMCNP(tallyDataList(counterDijColumns).name);
        disp(['Reading out tally data took: ', num2str(toc), ' seconds.'])
        disp('*****')

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

        % %% Calculate RMF parameters
        % if isfield(pln.propOpt,'bioOptimization')
        %     if strcmp(pln.propOpt.bioOptimization,'RBExSecPartDose_MCDS_RMFmodel')
        %         [dij, cst] = matRad_getRMFmodelParameters4neutrons(dij, counterDijColumns, doseMatrixBixel, pln, ct, cst);
        %     else
        %         %                     error('Biological optimization mode unknown for neutrons.')
        %     end
        % end


    end
end