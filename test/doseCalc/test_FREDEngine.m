function test_suite = test_FREDEngine

test_functions=localfunctions();

initTestSuite;

function test_constructFREDEngine
    radModes = DoseEngines.matRad_ParticleFREDEngine.possibleRadiationModes;
    for i = 1:numel(radModes)
        plnDummy = struct('radiationMode',radModes{i},'machine','Generic','propDoseCalc',struct('engine','FRED'));
        engine = DoseEngines.matRad_ParticleFREDEngine(plnDummy);
        assertTrue(isa(engine,'DoseEngines.matRad_ParticleFREDEngine'));
    end

function test_constructFailOnWrongRadMode
    plnDummy = struct('radiationMode','brachy','machine','HDR','propDoseCalc',struct('engine','FRED'));
    assertExceptionThrown(@()DoseEngines.matRad_ParticleFREDEngine(plnDummy));

function test_propertyAssignmentFromPln

    radModes = DoseEngines.matRad_ParticleFREDEngine.possibleRadiationModes;

    for i = 1:numel(radModes)
        pln = struct('radiationMode',radModes{i},'machine','Generic','propDoseCalc',struct('engine','FRED'));
        
        pln.propDoseCalc.HUclamping               = false;
        pln.propDoseCalc.HUtable                  = 'matRad_default_FredMaterialConverter';
        pln.propDoseCalc.exportCalculation        = true;
        pln.propDoseCalc.sourceModel              = 'gaussian';
        pln.propDoseCalc.useGPU                   = false;
        pln.propDoseCalc.roomMaterial             = 'Vacuum';
        pln.propDoseCalc.printOutput              = false;
        pln.propDoseCalc.numHistoriesDirect       = 42;
        pln.propDoseCalc.numHistoriesPerBeamlet   = 42;
        
        engine = DoseEngines.matRad_ParticleFREDEngine(pln);
        
        assertTrue(isa(engine,'DoseEngines.matRad_ParticleFREDEngine'));
        
        plnFields = fieldnames(pln.propDoseCalc);
        plnFields(strcmp([plnFields(:)], 'engine')) = [];
        
        for fieldIdx=1:numel(plnFields)
            assertTrue(isequal(engine.(plnFields{fieldIdx}), pln.propDoseCalc.(plnFields{fieldIdx})));
        end
    end

function test_writeFiles

        matRad_cfg = MatRad_Config.instance();
        radModes = DoseEngines.matRad_ParticleFREDEngine.possibleRadiationModes;

        load([radModes{1} '_testData.mat']);
        pln.radiationMode = radModes{1};
        pln.machine = 'Generic';
        pln.propDoseCalc.engine = 'FRED';
        pln.propDoseCalc.externalCalculation = 'write';

        w = ones(sum([stf(:).totalNumOfBixels]),1);

        resultGUI = matRad_calcDoseForward(ct,cst,stf,pln,w);

        fredMainDir   = [matRad_cfg.primaryUserFolder filesep 'FRED'];
        runFolder     = [fredMainDir filesep 'MCrun'];
        inputFolder   = [runFolder filesep 'inp'];
        planFolder    = [inputFolder filesep 'plan'];
        regionsFolder = [inputFolder filesep 'regions'];

        
        assertTrue(isfolder(fredMainDir));
        assertTrue(isfolder(runFolder));
        assertTrue(isfolder(inputFolder));
        assertTrue(isfolder(planFolder));
        assertTrue(isfolder(regionsFolder));        
        
        assertTrue(isfile([planFolder filesep 'plan.inp']));
        assertTrue(isfile([planFolder filesep 'planDelivery.inp']));

        assertTrue(isfile([regionsFolder filesep 'CTpatient.raw']));
        assertTrue(isfile([regionsFolder filesep 'CTpatient.mhd']));
        assertTrue(isfile([regionsFolder filesep 'regions.inp']));
        
        assertTrue(isfile([runFolder filesep 'fred.inp']));

function test_loadDij

        matRad_cfg = MatRad_Config.instance();
        load(['protons_testData.mat']);
        pln.machine = 'Generic';
        pln.propDoseCalc.engine = 'FRED';
        pln.propDoseCalc.useGPU = true;
        pln.propDoseCalc.externalCalculation = fullfile(matRad_cfg.matRadRoot, 'test', 'testData', 'FRED_data');
        
        % Test dij-load
        dijFredLoad         = matRad_calcDoseInfluence(ct,cst,stf,pln);
        
        % Test forward calculation cube load
        w = ones(sum([stf(:).totalNumOfBixels]),1);
        forwardDoseFredLoad = matRad_calcDoseForward(ct,cst,stf,pln,w);

        resultGUI = matRad_calcCubes(w, dijFredLoad, 1);

        nBixels = sum([stf(:).totalNumOfBixels]);
        nVoxles = prod(ct.cubeDim);

        % Assert basic parameters
        assertTrue(isequal(dijFredLoad.externalCalculationLodPath, fullfile(pln.propDoseCalc.externalCalculation, 'MCrun', 'out', 'scoreij', 'Phantom.Dose.bin')));
        assertTrue(isequal(size(dijFredLoad.physicalDose{1}),[nVoxles, nBixels]));
        assertTrue(isequal(size(forwardDoseFredLoad.physicalDose), size(resultGUI.physicalDose)));


function test_bioCalculation

        matRad_cfg = MatRad_Config.instance();
        load(['protons_testData.mat']);
        pln.machine = 'Generic';
        pln.propDoseCalc.bioModel = matRad_bioModel('protons', 'MCN');
        
        pln.propDoseCalc.engine = 'FRED';
        pln.propDoseCalc.useGPU = true;
        pln.propDoseCalc.externalCalculation = fullfile(matRad_cfg.matRadRoot, 'test', 'testData', 'FRED_data');
        
        % Test dij-load
        dijFredLoad         = matRad_calcDoseInfluence(ct,cst,stf,pln);
        
        % Test forward calculation cube load
        w = ones(sum([stf(:).totalNumOfBixels]),1);
        forwardDoseFredLoad = matRad_calcDoseForward(ct,cst,stf,pln,w);

        % Assert basic parameters
        assertTrue(all(cellfun(@(x) isfield(dijFredLoad, x), {'physicalDose', 'mLETd', 'mAlphaDose', 'mSqrtBetaDose'})));
        assertTrue(all(cellfun(@(x) isfield(forwardDoseFredLoad, x), {'physicalDose', 'LET', 'alpha', 'beta', 'effect'})));


function test_additionalParameters
        
        matRad_cfg = MatRad_Config.instance();
        load(['protons_testData.mat']);
        pln.machine = 'Generic';

        pln.propDoseCalc.HUtable = 'internal';
        pln.propDoseCalc.sourceModel = 'emittance';
        pln.propDoseCalc.HUclamping = true;

        engine = DoseEngines.matRad_ParticleFREDEngine(pln);

        assertTrue(all(cellfun(@(x,y) isequal(engine.(x), y), {'sourceModel','HUtable', 'HUclamping'}, {'emittance','internal', true})));