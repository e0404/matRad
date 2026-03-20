function test_suite = test_topasMCEngine

test_functions = localfunctions();

initTestSuite;

function test_loadMachine
radModes = DoseEngines.matRad_TopasMCEngine.possibleRadiationModes;
for i = 1:numel(radModes)
    machineName = 'Generic';
    machine = DoseEngines.matRad_TopasMCEngine.loadMachine(radModes{i}, machineName);
    assertTrue(isstruct(machine));
end
assertExceptionThrown(@() DoseEngines.matRad_TopasMCEngine.loadMachine('grbl', 'grbl'), 'matRad:Error');

function test_getEngineFromPlnByName
radModes = DoseEngines.matRad_TopasMCEngine.possibleRadiationModes;
for i = 1:numel(radModes)
    machineName = 'Generic';
    plnDummy = struct('radiationMode', radModes{i}, 'machine', machineName, 'propDoseCalc', struct('engine', 'TOPAS'));
    engine = DoseEngines.matRad_TopasMCEngine.getEngineFromPln(plnDummy);
    assertTrue(isa(engine, 'DoseEngines.matRad_TopasMCEngine'));
end

function test_topasMCdoseCalcBasic
% test if all the necessary output files are written for a couple of cases.
% i am not using the default number of histories for testing her, instead 1e6.
% Because the files are just written and not simulated so we dont care about simulation time.
% To few histories may result in weird behavior in the topas interface, i.e if a beam
% receives no histories because there are not enough to be distributed across the spots,
% it causes an error
radModes = DoseEngines.matRad_TopasMCEngine.possibleRadiationModes;
matRad_cfg = MatRad_Config.instance();

if moxunit_util_platform_is_octave
    confirm_recursive_rmdir(false, 'local');
end

for i = 1:numel(radModes)
    load([radModes{i} '_testData.mat']);
    pln.bioModel = matRad_bioModel(radModes{i}, 'none');

    w = ones(1, sum([stf(:).totalNumOfBixels]));

    if strcmp(radModes{i}, 'photons')
        pln.propOpt.runDAO = 1;
        dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
        resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels, 1), dij);
        resultGUI.wUnsequenced = ones(dij.totalNumOfBixels, 1);
        pln.propSeq.sequencer = 'siochi';
        pln.propSeq.sequencingLevel = 5;
        sequencer = matRad_SequencerBase.getSequencerFromPln(pln);
        resultGUI = matRad_sequencing(resultGUI, stf, pln);
        [pln, stf] = sequencer.aperture2collimation(pln, stf, resultGUI.sequencing);
        w = resultGUI.w;
        pln.propDoseCalc.beamProfile =  'phasespace';
    end

    pln.propDoseCalc.engine = 'TOPAS';
    pln.propDoseCalc.externalCalculation = 'write';
    pln.propDoseCalc.numHistoriesDirect = 1e6;
    resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, w);

    folderName = [matRad_cfg.primaryUserFolder filesep 'TOPAS' filesep];
    folderName = [folderName stf(1).radiationMode, '_', stf(1).machine, '_', datestr(now, 'dd-mm-yy')];
    % check of outputfolder exists
    assertTrue(isfolder(folderName));
    % check if file in folder existi
    assertTrue(isfile([folderName filesep 'matRad_cube.dat']));
    assertTrue(isfile([folderName filesep 'matRad_cube.txt']));
    assertTrue(isfile([folderName filesep 'MCparam.mat']));
    for j = 1:numel(stf)
        assertTrue(isfile([folderName filesep 'beamSetup_matRad_plan_field' num2str(j) '.txt']));
        assertTrue(isfile([folderName filesep 'matRad_plan_field' num2str(j) '_run1.txt']));
    end
    rmdir(folderName, 's'); % clean up
end

function test_topasMCdoseCalcBasicRBE
% test if all the necessary output files are written for a couple of cases.
% i am not using the default number of histories for testing her, instead 1e6.
% Because the files are just written and not simulated so we dont care about simulation time.
% To few histories may result in weird behavior in the topas interface, i.e if a beam
% receives no histories because there are not enough to be distributed across the spots,
% it causes an error
radModes = DoseEngines.matRad_TopasMCEngine.possibleRadiationModes;
matRad_cfg = MatRad_Config.instance();

if moxunit_util_platform_is_octave
    confirm_recursive_rmdir(false, 'local');
end

for i = 1:numel(radModes)
    switch radModes{i}
        case  'protons'
            RBEmodel = {'mcn', 'wed'};
        case {'helium', 'carbon'}
            RBEmodel = {'libamtrack', 'lem'};
        otherwise
            continue
    end
    matRad_cfg = MatRad_Config.instance();
    load([radModes{i} '_testData.mat']);

    pln.propDoseCalc.engine = 'TOPAS';
    pln.propDoseCalc.externalCalculation = 'write';
    pln.propDoseCalc.numHistoriesDirect = 1e6;
    pln.propDoseCalc.scorer.RBE = true;
    pln.propDoseCalc.scorer.RBE_model = RBEmodel;
    pln.bioModel = matRad_bioModel(radModes{i}, 'none');
    resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, ones(1, sum([stf(:).totalNumOfBixels])));

    folderName = [matRad_cfg.primaryUserFolder filesep 'TOPAS' filesep];
    folderName = [folderName stf(1).radiationMode, '_', stf(1).machine, '_', datestr(now, 'dd-mm-yy')];
    % check of outputfolder exists
    assertTrue(isfolder(folderName));
    % check if file in folder existi
    assertTrue(isfile([folderName filesep 'matRad_cube.dat']));
    assertTrue(isfile([folderName filesep 'matRad_cube.txt']));
    assertTrue(isfile([folderName filesep 'MCparam.mat']));
    for j = 1:numel(stf)
        assertTrue(isfile([folderName filesep 'beamSetup_matRad_plan_field' num2str(j) '.txt']));
        assertTrue(isfile([folderName filesep 'matRad_plan_field' num2str(j) '_run1.txt']));
    end
    rmdir(folderName, 's'); % clean up
end

function test_topasMCdoseCalcMultRuns
numOfRuns = 5;
radModes = DoseEngines.matRad_TopasMCEngine.possibleRadiationModes;
matRad_cfg = MatRad_Config.instance();

if moxunit_util_platform_is_octave
    confirm_recursive_rmdir(false, 'local');
end

for i = 1:numel(radModes)
    load([radModes{i} '_testData.mat']);
    pln.bioModel = matRad_bioModel(radModes{i}, 'none');
    w = ones(1, sum([stf(:).totalNumOfBixels]));

    if strcmp(radModes{i}, 'photons')
        dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
        resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels, 1), dij);
        resultGUI.wUnsequenced = ones(dij.totalNumOfBixels, 1);
        pln.propSeq.sequencer = 'siochi';
        pln.propSeq.sequencingLevel = 5;
        sequencer = matRad_SequencerBase.getSequencerFromPln(pln);
        resultGUI = matRad_sequencing(resultGUI, stf, pln);
        [pln, stf] = sequencer.aperture2collimation(pln, stf, resultGUI.sequencing);
        pln.propDoseCalc.beamProfile =  'phasespace';
        w = resultGUI.w;
    end

    pln.propDoseCalc.engine = 'TOPAS';
    pln.propDoseCalc.externalCalculation = 'write';
    pln.propDoseCalc.numHistoriesDirect = 1e6;
    pln.propDoseCalc.numOfRuns = numOfRuns;

    resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, w);

    folderName = [matRad_cfg.primaryUserFolder filesep 'TOPAS' filesep];
    folderName = [folderName stf(1).radiationMode, '_', stf(1).machine, '_', datestr(now, 'dd-mm-yy')];
    % check of outputfolder exists
    assertTrue(isfolder(folderName));
    % check if file in folder existi
    assertTrue(isfile([folderName filesep 'matRad_cube.dat']));
    assertTrue(isfile([folderName filesep 'matRad_cube.txt']));
    assertTrue(isfile([folderName filesep 'MCparam.mat']));
    for j = 1:numel(stf)
        assertTrue(isfile([folderName filesep 'beamSetup_matRad_plan_field' num2str(j) '.txt']));
        for k = 1:numOfRuns
            assertTrue(isfile([folderName filesep 'matRad_plan_field' num2str(j) '_run' num2str(k) '.txt']));
        end
    end

    rmdir(folderName, 's'); % clean up
end

function test_topasMCdoseCalc4D
numOfPhases = 5;
radModes = DoseEngines.matRad_TopasMCEngine.possibleRadiationModes;
matRad_cfg = MatRad_Config.instance();

if moxunit_util_platform_is_octave
    confirm_recursive_rmdir(false, 'local');
end

% physical Dose
for i = 1:numel(radModes)
    if ~ismember(radModes{i}, {'photons', 'VHEE'})
        load([radModes{i} '_testData.mat']);
        [ct, cst] = matRad_addMovement(ct, cst, 5, numOfPhases, [0 3 0], 'dvfType', 'pull');
        pln.bioModel = matRad_bioModel(radModes{i}, 'none');
        resultGUI.w = ones(1, sum([stf(:).totalNumOfBixels]))';
        resultGUI = matRad_sequencing(resultGUI, stf, pln);
        sequencer = matRad_ParticleSequencer();
        resultGUI.sequencing = sequencer.makePhaseMatrix(resultGUI.sequencing, ct.numOfCtScen, ct.motionPeriod);
        pln.propDoseCalc.engine = 'TOPAS';
        pln.propDoseCalc.externalCalculation = 'write';
        pln.propDoseCalc.calc4DInterplay = true;
        pln.propDoseCalc.calcTimeSequence = resultGUI.sequencing;
        pln.propDoseCalc.numHistoriesDirect = 1e6;
        resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, resultGUI.w);

        folderName = [matRad_cfg.primaryUserFolder filesep 'TOPAS' filesep];
        folderName = [folderName stf(1).radiationMode, '_', stf(1).machine, '_', datestr(now, 'dd-mm-yy')];
        % check of outputfolder exists
        assertTrue(isfolder(folderName));
        % check if file in folder existi
        assertTrue(isfile([folderName filesep 'MCparam.mat']));
        for j = 1:numel(stf)
            assertTrue(isfile([folderName filesep 'beamSetup_matRad_plan_field' num2str(j) '.txt']));
            assertTrue(isfile([folderName filesep 'matRad_plan_field' num2str(j) '_run1.txt']));
            assertTrue(isfile([folderName filesep 'matRad_cube_field' num2str(j) '.txt']));
            for k = 1:numOfPhases
                assertTrue(isfile([folderName filesep 'matRad_cube' num2str(k) '.dat']));
            end
        end
        rmdir(folderName, 's'); % clean up
    end
end
% RBExDose
for i = 1:numel(radModes)
    switch radModes{i}
        case  'protons'
            RBEmodel = {'mcn', 'wed'};
        case {'helium', 'carbon'}
            RBEmodel = {'libamtrack', 'lem'};
        otherwise
            continue
    end

    load([radModes{i} '_testData.mat']);
    [ct, cst] = matRad_addMovement(ct, cst, 5, numOfPhases, [0 3 0], 'dvfType', 'pull');
    pln.bioModel = matRad_bioModel(radModes{i}, 'none');
    resultGUI.w = ones(1, sum([stf(:).totalNumOfBixels]))';
    resultGUI = matRad_sequencing(resultGUI, stf, pln);
    sequencer = matRad_ParticleSequencer();
    resultGUI.sequencing = sequencer.makePhaseMatrix(resultGUI.sequencing, ct.numOfCtScen, ct.motionPeriod);
    pln.propDoseCalc.engine = 'TOPAS';
    pln.propDoseCalc.externalCalculation = 'write';
    pln.propDoseCalc.calc4DInterplay = true;
    pln.propDoseCalc.calcTimeSequence = resultGUI.sequencing;
    pln.propDoseCalc.numHistoriesDirect = 1e6;
    pln.propDoseCalc.scorer.RBE = true;
    pln.propDoseCalc.scorer.RBE_model = RBEmodel;
    resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, resultGUI.w);

    folderName = [matRad_cfg.primaryUserFolder filesep 'TOPAS' filesep];
    folderName = [folderName stf(1).radiationMode, '_', stf(1).machine, '_', datestr(now, 'dd-mm-yy')];
    % check of outputfolder exists
    assertTrue(isfolder(folderName));
    % check if file in folder existi
    assertTrue(isfile([folderName filesep 'MCparam.mat']));
    for j = 1:numel(stf)
        assertTrue(isfile([folderName filesep 'beamSetup_matRad_plan_field' num2str(j) '.txt']));
        assertTrue(isfile([folderName filesep 'matRad_plan_field' num2str(j) '_run1.txt']));
        assertTrue(isfile([folderName filesep 'matRad_cube_field' num2str(j) '.txt']));
        for k = 1:numOfPhases
            assertTrue(isfile([folderName filesep 'matRad_cube' num2str(k) '.dat']));
        end
    end
    rmdir(folderName, 's'); % clean up
end

function test_topasMultiAlphaBeta
radModes = {'protons'};
matRad_cfg = MatRad_Config.instance();

if moxunit_util_platform_is_octave
    confirm_recursive_rmdir(false, 'local');
end

for i = 1:numel(radModes)
    switch radModes{i}
        case  'protons'
            RBEmodel = {'mcn', 'wed'};
        case {'helium', 'carbon'}
            RBEmodel = {'libamtrack', 'lem'};
        otherwise
            continue
    end
    matRad_cfg = MatRad_Config.instance();
    load([radModes{i} '_testData.mat']);

    % Extract purely random alpha/beta couples
    for idx = 1:size(cst, 1)
        alphaList(1, idx) =  0.1 + (0.8 - 0.1) .* rand;
        betaList(1, idx) = 0.001 + (0.1 - 0.001) .* rand;
    end
    % alpha1 = 0.1 + (0.8-0.1).*rand;
    % alpha2 = 0.1 + (0.8-0.1).*rand;
    % beta1  = 0.001 + (0.1-0.001).*rand;
    % beta2  = 0.001 + (0.1-0.001).*rand;
    rbeIdx = 1 + (length(RBEmodel) - 1) * round(rand);

    for idx = 1:size(cst, 1)
        cst{idx, 5}.alphaX = alphaList(idx);
        cst{idx, 5}.betaX = betaList(idx);
        % cst{1,5}.alphaX = alpha1;
        % cst{1,5}.betaX = beta1;
        % cst{2,5}.alphaX = alpha2;
        % cst{2,5}.betaX = beta2;
    end

    pln.propDoseCalc.bioParameters.AlphaX = alphaList; % [alpha1, alpha2];
    pln.propDoseCalc.bioParameters.BetaX = betaList; % [beta1, beta2];

    pln.propDoseCalc.engine                 = 'TOPAS';
    pln.propDoseCalc.externalCalculation    = 'write';
    pln.propDoseCalc.numHistoriesDirect     = 1e6;
    pln.propDoseCalc.numOfRuns              = 1;
    pln.propDoseCalc.scorer.RBE             = true;
    pln.propDoseCalc.scorer.RBE_model       = {RBEmodel{rbeIdx}};
    % pln.bioModel = matRad_bioModel(radModes{i},'none');
    pln.bioModel = upper(RBEmodel{rbeIdx});
    pln.propDoseCalc.scorer.reportQuantity  = {'Sum'};
    resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, ones(1, sum([stf(:).totalNumOfBixels])));
    % matRad_calcDoseForward(ct,cst,stf,pln,resultGUI.w);

    folderName = [matRad_cfg.primaryUserFolder filesep 'TOPAS' filesep];
    folderName = [folderName stf(1).radiationMode, '_', stf(1).machine, '_', datestr(now, 'dd-mm-yy')];
    % check of outputfolder exists
    assertTrue(isfolder(folderName));
    % check if file in folder existi
    assertTrue(isfile([folderName filesep 'matRad_cube.dat']));
    assertTrue(isfile([folderName filesep 'matRad_cube.txt']));
    assertTrue(isfile([folderName filesep 'MCparam.mat']));
    for j = 1:pln.propStf.numOfBeams
        assertTrue(isfile([folderName filesep 'beamSetup_matRad_plan_field' num2str(j) '.txt']));
        assertTrue(isfile([folderName filesep 'matRad_plan_field' num2str(j) '_run1.txt']));

        for runIdx = 1:pln.propDoseCalc.numOfRuns
            pathFile = [folderName filesep 'matRad_plan_field' num2str(j) '_run' num2str(runIdx) '.txt'];
            fileText = fileread(pathFile);
            matches = unique(regexp(fileText, 'CellType_\d+/', 'match'));
            assertTrue(length(matches) == length(alphaList) && length(matches) == length(betaList));

            lines = strsplit(fileText, '\n');
            idx = strncmp(lines, 'd:Sc/AlphaX', length('d:Sc/AlphaX'));
            matchingLines = lines(idx);
            for cellIdx = 1:length(matchingLines)
                tmp = strsplit(matchingLines{cellIdx}, '_');
                tmp = strsplit(tmp{end}, ' ');
                assertTrue(abs(alphaList(str2double(tmp{1})) - str2double(tmp{3})) <= 1e-4);
            end

            lines = strsplit(fileText, '\n');
            idx = strncmp(lines, 'd:Sc/BetaX', length('d:Sc/BetaX'));
            matchingLines = lines(idx);
            for cellIdx = 1:length(matchingLines)
                tmp = strsplit(matchingLines{cellIdx}, '_');
                tmp = strsplit(tmp{end}, ' ');
                assertTrue(abs(betaList(str2double(tmp{1})) - str2double(tmp{3})) < 1e-4);
            end

        end

    end
end
