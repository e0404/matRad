function test_suite = test_Analytical

test_functions=localfunctions();

initTestSuite;

function test_getAnalyticalEngineFromPln
    % Single gaussian lateral model
    testData.pln = struct('radiationMode','protons','machine','Generic');
    testData.pln.propDoseCalc.engine = 'AnalyticalPB';
    engine = DoseEngines.matRad_ParticleAnalyticalBortfeldEngine.getEngineFromPln(testData.pln);
    assertTrue(isa(engine,'DoseEngines.matRad_ParticleAnalyticalBortfeldEngine'));

function test_loadMachineForAnalytical
    possibleRadModes = DoseEngines.matRad_ParticleAnalyticalBortfeldEngine.possibleRadiationModes;
    for i = 1:numel(possibleRadModes)
        machine = DoseEngines.matRad_ParticleAnalyticalBortfeldEngine.loadMachine(possibleRadModes{i},'Generic');
        assertTrue(isstruct(machine));
        assertTrue(isfield(machine, 'meta'));
        assertTrue(isfield(machine.meta, 'radiationMode'));
        assertTrue(strcmp(machine.meta.radiationMode, 'protons'));
    end

function test_calcDoseAnalytical
    testData = load('protons_testData.mat');
    assertTrue(DoseEngines.matRad_ParticleAnalyticalBortfeldEngine.isAvailable(testData.pln))
    testData.pln.propDoseCalc.engine = 'AnalyticalPB';
    resultGUI = matRad_calcDoseForward(testData.ct, testData.cst, testData.stf, testData.pln, ones(sum([testData.stf(:).totalNumOfBixels]),1));

    assertTrue(isfield(resultGUI, 'physicalDose'));
    assertTrue(isfield(resultGUI, 'w'));
    assertTrue(isequal(testData.ct.cubeDim, size(resultGUI.physicalDose)));

function test_nonSupportedSettings
    % Radiation mode other than protons not implemented 
    testData = load('carbon_testData.mat');
    testData.pln.propDoseCalc.engine = 'AnalyticalPB';    
    assertFalse(DoseEngines.matRad_ParticleAnalyticalBortfeldEngine.isAvailable(testData.pln));
   
    % Biological models, LET, other lateral models not implemented
    testData = load('protons_testData.mat'); 
    testData.pln.propDoseCalc.engine = 'AnalyticalPB';
    testData.pln.propDoseCalc.calcLET = true;
    testData.pln.propDoseCalc.calcBioDose = true;
    testData.pln.propDoseCalc.lateralModel = 'double';
    engine = DoseEngines.matRad_ParticleAnalyticalBortfeldEngine.getEngineFromPln(testData.pln);
    assertTrue(isa(engine,'DoseEngines.matRad_ParticleAnalyticalBortfeldEngine'));

    resultGUI = matRad_calcDoseForward(testData.ct, testData.cst, testData.stf, testData.pln, ones(sum([testData.stf(:).totalNumOfBixels]),1));
    assertTrue(~engine.calcLET)
    %assertTrue(~engine.calcBioDose) % Access protected property

    % Invalid machine without radiation mode field
    testData.pln.machine = 'Empty';
    testData.pln.propDoseCalc.engine = 'AnalyticalPB';    
    assertExceptionThrown(@() DoseEngines.matRad_ParticleAnalyticalBortfeldEngine.isAvailable(testData.pln));
    assertFalse(DoseEngines.matRad_ParticleAnalyticalBortfeldEngine.isAvailable(testData.pln,[]));



