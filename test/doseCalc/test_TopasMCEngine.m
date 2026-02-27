function test_suite = test_TopasMCEngine

test_functions=localfunctions();

initTestSuite;

function test_loadMachine
    radModes = DoseEngines.matRad_TopasMCEngine.possibleRadiationModes;
    for i = 1:numel(radModes)
        machineName = 'Generic';
        machine = DoseEngines.matRad_TopasMCEngine.loadMachine(radModes{i},machineName);
        assertTrue(isstruct(machine));
    end
    assertExceptionThrown(@() DoseEngines.matRad_TopasMCEngine.loadMachine('grbl','grbl'),'matRad:Error')

function test_getEngineFromPlnByName
    radModes = DoseEngines.matRad_TopasMCEngine.possibleRadiationModes;
    for i = 1:numel(radModes)
        machineName = 'Generic';
        plnDummy = struct('radiationMode',radModes{i},'machine',machineName,'propDoseCalc',struct('engine','TOPAS'));
         engine = DoseEngines.matRad_TopasMCEngine.getEngineFromPln(plnDummy);
        assertTrue(isa(engine,'DoseEngines.matRad_TopasMCEngine'));
    end

function test_TopasMCdoseCalcBasic
% test if all the necessary output files are written vor a couple of cases.
% i am not using the default number of histories for testing her, insted 1e6.
% Because the files are just writen and not simulated so we dont care about simulation time. 
% To few histories may result in wierd behavior in the topas interface, i.e if a beam
% recieves no histories because there are not enough to be distributed accros the spots,
% it causes an error
radModes = DoseEngines.matRad_TopasMCEngine.possibleRadiationModes;
matRad_cfg = MatRad_Config.instance();

if moxunit_util_platform_is_octave
  confirm_recursive_rmdir(false,'local');
end

for i = 1:numel(radModes)
    load([radModes{i} '_testData.mat']);
    pln.bioModel = matRad_bioModel(radModes{i},'none');
    
    w = ones(1,sum([stf(:).totalNumOfBixels]));
    
    if strcmp(radModes{i},'photons')
        pln.propOpt.runSequencing =  1;
        pln.propOpt.runDAO = 1;
        dij = matRad_calcDoseInfluence(ct,cst,stf,pln);
        resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels,1),dij);
        resultGUI.wUnsequenced = ones(dij.totalNumOfBixels,1);
        resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,5);
        [pln,stf] = matRad_aperture2collimation(pln,stf,resultGUI.sequencing,resultGUI.apertureInfo);
        w = resultGUI.w;
        pln.propDoseCalc.beamProfile =  'phasespace';
    end

    pln.propDoseCalc.engine = 'TOPAS';
    pln.propDoseCalc.externalCalculation = 'write';
    pln.propDoseCalc.numHistoriesDirect = 1e6;
    resultGUI = matRad_calcDoseForward(ct,cst,stf,pln, w);
    
    folderName = [matRad_cfg.primaryUserFolder filesep 'TOPAS' filesep];
    folderName = [folderName stf(1).radiationMode,'_',stf(1).machine,'_',datestr(now, 'dd-mm-yy')];
    %check of outputfolder exists
    assertTrue(isfolder(folderName));
    %check if file in folder existi
    assertTrue(isfile([folderName filesep 'matRad_cube.dat']));
    assertTrue(isfile([folderName filesep 'matRad_cube.txt']));
    assertTrue(isfile([folderName filesep 'MCparam.mat']));
    for j = 1:pln.propStf.numOfBeams
        assertTrue(isfile([folderName filesep 'beamSetup_matRad_plan_field' num2str(j) '.txt']));
        assertTrue(isfile([folderName filesep 'matRad_plan_field' num2str(j) '_run1.txt']));
    end
    rmdir(folderName,'s'); %clean up
end


function test_TopasMCdoseCalcBasicRBE
% test if all the necessary output files are written vor a couple of cases.
% i am not using the default number of histories for testing her, insted 1e6.
% Because the files are just writen and not simulated so we dont care about simulation time. 
% To few histories may result in wierd behavior in the topas interface, i.e if a beam
% recieves no histories because there are not enough to be distributed accros the spots,
% it causes an error
radModes = DoseEngines.matRad_TopasMCEngine.possibleRadiationModes;
matRad_cfg = MatRad_Config.instance();

if moxunit_util_platform_is_octave
  confirm_recursive_rmdir(false,'local');
end

for i = 1:numel(radModes)
    switch radModes{i}
        case  'protons'
            RBEmodel = {'mcn', 'wed'};
        case {'helium', 'carbon'} 
            RBEmodel ={'libamtrack','lem'};
        otherwise
            continue;
    end
    matRad_cfg = MatRad_Config.instance();
    load([radModes{i} '_testData.mat']);
    
    pln.propDoseCalc.engine = 'TOPAS';
    pln.propDoseCalc.externalCalculation = 'write';
    pln.propDoseCalc.numHistoriesDirect = 1e6;
    pln.propDoseCalc.scorer.RBE = true;
    pln.propDoseCalc.scorer.RBE_model = RBEmodel;
    pln.bioModel = matRad_bioModel(radModes{i},'none');
    resultGUI = matRad_calcDoseForward(ct,cst,stf,pln, ones(1,sum([stf(:).totalNumOfBixels])));

    folderName = [matRad_cfg.primaryUserFolder filesep 'TOPAS' filesep];
    folderName = [folderName stf(1).radiationMode,'_',stf(1).machine,'_',datestr(now, 'dd-mm-yy')];
    %check of outputfolder exists
    assertTrue(isfolder(folderName));
    %check if file in folder existi
    assertTrue(isfile([folderName filesep 'matRad_cube.dat']));
    assertTrue(isfile([folderName filesep 'matRad_cube.txt']));
    assertTrue(isfile([folderName filesep 'MCparam.mat']));
    for j = 1:pln.propStf.numOfBeams
        assertTrue(isfile([folderName filesep 'beamSetup_matRad_plan_field' num2str(j) '.txt']));
        assertTrue(isfile([folderName filesep 'matRad_plan_field' num2str(j) '_run1.txt']));
    end
    rmdir(folderName,'s'); %clean up
end

function test_TopasMCdoseCalcMultRuns
numOfRuns = 5;
radModes = DoseEngines.matRad_TopasMCEngine.possibleRadiationModes;
matRad_cfg = MatRad_Config.instance();

if moxunit_util_platform_is_octave
  confirm_recursive_rmdir(false,'local');
end

for i = 1:numel(radModes)
    load([radModes{i} '_testData.mat']);
    pln.bioModel = matRad_bioModel(radModes{i},'none');
    w = ones(1,sum([stf(:).totalNumOfBixels]));
    
    if strcmp(radModes{i},'photons')
        pln.propOpt.runSequencing =  1;
        pln.propOpt.runDAO = 1;
        dij = matRad_calcDoseInfluence(ct,cst,stf,pln);
        resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels,1),dij);
        resultGUI.wUnsequenced = ones(dij.totalNumOfBixels,1);
        resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,5);
        [pln,stf] = matRad_aperture2collimation(pln,stf,resultGUI.sequencing,resultGUI.apertureInfo);
        pln.propDoseCalc.beamProfile =  'phasespace';
        w = resultGUI.w;
    end
    
    pln.propDoseCalc.engine = 'TOPAS';
    pln.propDoseCalc.externalCalculation = 'write';
    pln.propDoseCalc.numHistoriesDirect = 1e6;
    pln.propDoseCalc.numOfRuns = numOfRuns;

    resultGUI = matRad_calcDoseForward(ct,cst,stf,pln, w);

    folderName = [matRad_cfg.primaryUserFolder filesep 'TOPAS' filesep];
    folderName = [folderName stf(1).radiationMode,'_',stf(1).machine,'_',datestr(now, 'dd-mm-yy')];
    %check of outputfolder exists
    assertTrue(isfolder(folderName));
    %check if file in folder existi
    assertTrue(isfile([folderName filesep 'matRad_cube.dat']));
    assertTrue(isfile([folderName filesep 'matRad_cube.txt']));
    assertTrue(isfile([folderName filesep 'MCparam.mat']));
    for j = 1:pln.propStf.numOfBeams
        assertTrue(isfile([folderName filesep 'beamSetup_matRad_plan_field' num2str(j) '.txt']));
        for k = 1:numOfRuns
            assertTrue(isfile([folderName filesep 'matRad_plan_field' num2str(j) '_run' num2str(k) '.txt']));
        end
    end

    rmdir(folderName,'s'); %clean up
end


function test_TopasMCdoseCalc4D
numOfPhases = 5;
radModes = DoseEngines.matRad_TopasMCEngine.possibleRadiationModes;
matRad_cfg = MatRad_Config.instance();

if moxunit_util_platform_is_octave
  confirm_recursive_rmdir(false,'local');
end

% physical Dose
for i = 1:numel(radModes)
    if ~strcmp(radModes{i},'photons')
        load([radModes{i} '_testData.mat']);
        [ct,cst] = matRad_addMovement(ct, cst,5, numOfPhases,[0 3 0],'dvfType','pull');
        pln.bioModel = matRad_bioModel(radModes{i},'none');
        resultGUI.w = ones(1,sum([stf(:).totalNumOfBixels]))';
        timeSequence = matRad_makeBixelTimeSeq(stf, resultGUI);
        timeSequence = matRad_makePhaseMatrix(timeSequence, ct.numOfCtScen, ct.motionPeriod, 'linear');
        pln.propDoseCalc.engine = 'TOPAS';
        pln.propDoseCalc.externalCalculation = 'write';
        pln.propDoseCalc.calc4DInterplay = true;
        pln.propDoseCalc.calcTimeSequence = timeSequence;
        pln.propDoseCalc.numHistoriesDirect = 1e6;
        resultGUI = matRad_calcDoseForward(ct,cst,stf,pln, resultGUI.w);

        folderName = [matRad_cfg.primaryUserFolder filesep 'TOPAS' filesep];
        folderName = [folderName stf(1).radiationMode,'_',stf(1).machine,'_',datestr(now, 'dd-mm-yy')];
        %check of outputfolder exists
        assertTrue(isfolder(folderName));
        %check if file in folder existi
        assertTrue(isfile([folderName filesep 'MCparam.mat']));
        for j = 1:pln.propStf.numOfBeams
            assertTrue(isfile([folderName filesep 'beamSetup_matRad_plan_field' num2str(j) '.txt']));
            assertTrue(isfile([folderName filesep 'matRad_plan_field' num2str(j) '_run1.txt']));
            assertTrue(isfile([folderName filesep 'matRad_cube_field' num2str(j) '.txt']));
            for k = 1:numOfPhases
                assertTrue(isfile([folderName filesep 'matRad_cube' num2str(k) '.dat']));
            end
        end
        rmdir(folderName,'s'); %clean up
    end
end
%RBExDose
for i = 1:numel(radModes)
    switch radModes{i}
        case  'protons'
            RBEmodel = {'mcn', 'wed'};
        case {'helium', 'carbon'} 
            RBEmodel ={'libamtrack','lem'};
        otherwise
            continue;
    end

    load([radModes{i} '_testData.mat']);
    [ct,cst] = matRad_addMovement(ct, cst,5, numOfPhases,[0 3 0],'dvfType','pull');
    pln.bioModel = matRad_bioModel(radModes{i},'none');
    resultGUI.w = ones(1,sum([stf(:).totalNumOfBixels]))';
    timeSequence = matRad_makeBixelTimeSeq(stf, resultGUI);
    timeSequence = matRad_makePhaseMatrix(timeSequence, ct.numOfCtScen, ct.motionPeriod, 'linear');
    pln.propDoseCalc.engine = 'TOPAS';
    pln.propDoseCalc.externalCalculation = 'write';
    pln.propDoseCalc.calc4DInterplay = true;
    pln.propDoseCalc.calcTimeSequence = timeSequence;
    pln.propDoseCalc.numHistoriesDirect = 1e6;
    pln.propDoseCalc.scorer.RBE = true;
    pln.propDoseCalc.scorer.RBE_model = RBEmodel;
    resultGUI = matRad_calcDoseForward(ct,cst,stf,pln, resultGUI.w);

    folderName = [matRad_cfg.primaryUserFolder filesep 'TOPAS' filesep];
    folderName = [folderName stf(1).radiationMode,'_',stf(1).machine,'_',datestr(now, 'dd-mm-yy')];
    %check of outputfolder exists
    assertTrue(isfolder(folderName));
    %check if file in folder existi
    assertTrue(isfile([folderName filesep 'MCparam.mat']));
    for j = 1:pln.propStf.numOfBeams
        assertTrue(isfile([folderName filesep 'beamSetup_matRad_plan_field' num2str(j) '.txt']));
        assertTrue(isfile([folderName filesep 'matRad_plan_field' num2str(j) '_run1.txt']));
        assertTrue(isfile([folderName filesep 'matRad_cube_field' num2str(j) '.txt']));
        for k = 1:numOfPhases
            assertTrue(isfile([folderName filesep 'matRad_cube' num2str(k) '.dat']));
        end
    end
    rmdir(folderName,'s'); %clean up
end



