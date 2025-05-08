function dij = matRad_evaluateTallyMCNP(dij, cst, ct)
%% Read output from MCNP calculation and generate dij matix
%% Preparation and get list of mctal data
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

%% Resize cst (and forget about it)
cst = matRad_resizeCstToGrid(cst,ct.resolution.y:ct.resolution.y:ct.cubeDim(1)*ct.resolution.y,...
    ct.resolution.x:ct.resolution.x:ct.cubeDim(2)*ct.resolution.x,...
    ct.resolution.z:ct.resolution.z:ct.cubeDim(3)*ct.resolution.z,...
    dij.doseGrid.resolution.y:dij.doseGrid.resolution.y:dij.doseGrid.dimensions(1)*dij.doseGrid.resolution.y,...
    dij.doseGrid.resolution.x:dij.doseGrid.resolution.x:dij.doseGrid.dimensions(2)*dij.doseGrid.resolution.x,...
    dij.doseGrid.resolution.z:dij.doseGrid.resolution.z:dij.doseGrid.dimensions(3)*dij.doseGrid.resolution.z);

%% Allocate dij dose matrix 
dij.physicalDose{1} = spalloc(prod(dij.doseGrid.dimensions),dij.totalNumOfBixels,1);
for counterBixel=1:dij.totalNumOfBixels
    dij.doseMatrixBixel(counterBixel).physicalDose_relError{1} = spalloc(prod(dij.doseGrid.dimensions),1,1);
end

%% Read TMESH results
counterDijColumns = 0; % counter to go through the columns of dij

for counterBeam = 1:dij.numOfBeams
    for counterBixel = 1:dij.numOfRaysPerBeam(counterBeam)
        counterDijColumns = counterDijColumns+1;

        dij.bixelNum(counterDijColumns) = counterBixel;
        dij.rayNum(counterDijColumns) = counterBixel;
        dij.beamNum(counterDijColumns) = counterBeam;

        % Evaluate mesh tally and store neutron and photon dose
        matRad_cfg.dispInfo('*****\n')
        matRad_cfg.dispInfo('Evaluate Tally Data.\n')
        matRad_cfg.dispInfo('*****\n')
        tic;

        % Read TMESH results
        resultMCNP = matRad_readDataFromText_TMESHvBioOpti(dummyList(counterDijColumns).name, 'TMESH3', 2);
        resultMCNP = resultMCNP';
        doseMatrixBixel.physicalDose = zeros(dij.doseGrid.dimensions(2), dij.doseGrid.dimensions(1), dij.doseGrid.dimensions(3));   % Total dose
        doseMatrixBixel.physicalDose(1:end) = resultMCNP(:,1);


        doseMatrixBixel.physicalDose_relError = zeros(dij.doseGrid.dimensions(2), dij.doseGrid.dimensions(1), dij.doseGrid.dimensions(3));  % Relative error of total dose
        doseMatrixBixel.physicalDose_relError(1:end) = resultMCNP(:,2);

        clear resultMCNP

        doseMatrixBixel.physicalDose = permute(doseMatrixBixel.physicalDose, [2,1,3]);    % Permute matrix to match matRad coordinate system
        for cellCounter = 1:size(ct.doseGridCT.tissueBin,2)
            if cellCounter == 1
                doseMatrixBixel.physicalDose(ct.doseGridCT.tissueBin(cellCounter).linIndVol) = 0;   % Dose deposition in air is neglected
            else
                doseMatrixBixel.physicalDose(ct.doseGridCT.tissueBin(cellCounter).linIndVol) = ...
                    doseMatrixBixel.physicalDose(ct.doseGridCT.tissueBin(cellCounter).linIndVol)./mean(ct.doseGridCT.density{1,1}(ct.doseGridCT.tissueBin(cellCounter).linIndVol));   % MCNP output is in MeV/cm^3/source particle & ct.doseGridCT.density is given in g/cm^3
            end
            % mean(ct.doseGridCT.density{1,1}(ct.doseGridCT.tissueBin(cellCounter).linIndVol))
        end
        doseMatrixBixel.physicalDose = doseMatrixBixel.physicalDose*1.602177e-19*1e6*1e3; % Convert MeV/g to J/kg, output is now in Gy/source particle
        doseMatrixBixel.physicalDose = doseMatrixBixel.physicalDose./max(doseMatrixBixel.physicalDose, [], 'all');

        doseMatrixBixel.physicalDose_relError = permute(doseMatrixBixel.physicalDose_relError, [2,1,3]);

        % Bad TMESH tally statistics can lead to negative results
        if sum(doseMatrixBixel.physicalDose<0, 'all')
            matRad_cfg.dispWarning('*********************\n')
            matRad_cfg.dispWarning('Negative TMESH tally results detected. This is a hint for bad statistics!')
            matRad_cfg.dispWarning('Negative results set to zero.')
            matRad_cfg.dispWarning(['Minimum value: ', num2str(min(doseMatrixBixel.physicalDose, [], 'all')), ' Gy and maximum value: ', num2str(max(doseMatrixBixel.physicalDose, [], 'all')), ' Gy.'])
            matRad_cfg.dispWarning('*********************')

            doseMatrixBixel.physicalDose(doseMatrixBixel.physicalDose<0) = 0;
            doseMatrixBixel.physicalDose_relError(doseMatrixBixel.physicalDose<0) = 0;
        end

        matRad_cfg.dispInfo(['Reading out tally data took: ', num2str(toc), ' seconds.\n'])
        matRad_cfg.dispInfo('*****\n')

        dij.doseMatrixBixel(counterDijColumns).physicalDose_relError{1} = ...
            sparse(find(doseMatrixBixel.physicalDose_relError),1,doseMatrixBixel.physicalDose_relError(doseMatrixBixel.physicalDose_relError~=0),dij.doseGrid.numOfVoxels,1);

        % Store total dose from TMESH tally
        dij.physicalDose{1}(:,counterDijColumns) = ...
            sparse(find(doseMatrixBixel.physicalDose), 1, doseMatrixBixel.physicalDose(doseMatrixBixel.physicalDose~=0), dij.doseGrid.numOfVoxels,1);

        % Save relative error for each radiotherapy structure and bixel
        for counterRTStruct=1:size(cst,1)
            dij.relError_radioTherpyStruct(counterRTStruct,counterDijColumns).name = cst{counterRTStruct,2};
            dij.relError_radioTherpyStruct(counterRTStruct,counterDijColumns).bixelNumber = counterDijColumns;
            dij.relError_radioTherpyStruct(counterRTStruct,counterDijColumns).meanError = mean(doseMatrixBixel.physicalDose_relError(cst{counterRTStruct,4}{1,1}));
            dij.relError_radioTherpyStruct(counterRTStruct,counterDijColumns).maxError = max(doseMatrixBixel.physicalDose_relError(cst{counterRTStruct,4}{1,1}));
            dij.relError_radioTherpyStruct(counterRTStruct,counterDijColumns).medianError = median(doseMatrixBixel.physicalDose_relError(cst{counterRTStruct,4}{1,1}));
        end

        %% Calculate RMF parameters
        % if isfield(pln.propOpt,'bioOptimization')
        %     if strcmp(pln.propOpt.bioOptimization,'RBExSecPartDose_MCDS_RMFmodel')
        %         [dij, cst] = matRad_getRMFmodelParameters4neutrons(dij, counterDijColumns, doseMatrixBixel, pln, ct, cst);
        %     else
        %         %                     error('Biological optimization mode unknown for neutrons.')
        %     end
        % end


    end
end