function test_suite = test_HongPB

    test_functions=localfunctions();
    
    initTestSuite;
    
    function test_getHongPBEngineFromPln
        % Single gaussian lateral model
        testData.pln = struct('radiationMode','protons','machine','Generic');
        testData.pln.propDoseCalc.engine = 'HongPB';
        engine = DoseEngines.matRad_ParticleHongPencilBeamEngine.getEngineFromPln(testData.pln);
        assertTrue(isa(engine,'DoseEngines.matRad_ParticleHongPencilBeamEngine'));
        
    function test_loadMachineForHongPB
        possibleRadModes = DoseEngines.matRad_ParticleHongPencilBeamEngine.possibleRadiationModes;
        for i = 1:numel(possibleRadModes)
            machine = DoseEngines.matRad_ParticleHongPencilBeamEngine.loadMachine(possibleRadModes{i},'Generic');
            assertTrue(isstruct(machine));
            assertTrue(isfield(machine, 'meta'));
            assertTrue(isfield(machine.meta, 'radiationMode'));
            assertTrue(strcmp(machine.meta.radiationMode, possibleRadModes{i}));
        end
    
    function test_calcDoseHongPBprotons
        testData = load('protons_testData.mat');

        assertTrue(DoseEngines.matRad_ParticleHongPencilBeamEngine.isAvailable(testData.pln));
        
        testData.pln.propDoseCalc.engine = 'HongPB';
        testData.pln.propDoseCalc.dosimetricLateralCutOff = 0.995;
        testData.pln.propDoseCalc.geometricLateralCutOff = 50;
        resultGUI = matRad_calcDoseForward(testData.ct, testData.cst, testData.stf, testData.pln, ones(sum([testData.stf(:).totalNumOfBixels]),1));
    
        assertTrue(isequal(fieldnames(resultGUI),fieldnames(testData.resultGUI)));
        assertTrue(isequal(testData.ct.cubeDim, size(resultGUI.physicalDose)));
        assertElementsAlmostEqual(resultGUI.physicalDose,testData.resultGUI.physicalDose,'relative',1e-2,1e-2);

    function test_calcDoseHongPBhelium
        testData = load('helium_testData.mat');
        assertTrue(DoseEngines.matRad_ParticleHongPencilBeamEngine.isAvailable(testData.pln));

        testData.pln.propDoseCalc.engine = 'HongPB';
        testData.pln.propDoseCalc.dosimetricLateralCutOff = 0.995;
        testData.pln.propDoseCalc.geometricLateralCutOff = 50;
        resultGUI = matRad_calcDoseForward(testData.ct, testData.cst, testData.stf, testData.pln, ones(sum([testData.stf(:).totalNumOfBixels]),1));
    
        assertTrue(isequal(fieldnames(resultGUI),fieldnames(testData.resultGUI)));
        assertTrue(isequal(testData.ct.cubeDim, size(resultGUI.physicalDose)));
        assertElementsAlmostEqual(resultGUI.physicalDose,testData.resultGUI.physicalDose,'relative',1e-2,1e-2);

    function test_calcDoseHongPBcarbon
        testData = load('carbon_testData.mat');
        assertTrue(DoseEngines.matRad_ParticleHongPencilBeamEngine.isAvailable(testData.pln));

        testData.pln.propDoseCalc.engine = 'HongPB';
        testData.pln.propDoseCalc.dosimetricLateralCutOff = 0.995;
        testData.pln.propDoseCalc.geometricLateralCutOff = 50;
        resultGUI = matRad_calcDoseForward(testData.ct, testData.cst, testData.stf, testData.pln, ones(sum([testData.stf(:).totalNumOfBixels]),1));
    
        assertTrue(isequal(fieldnames(resultGUI),fieldnames(testData.resultGUI)));
        assertTrue(isequal(testData.ct.cubeDim, size(resultGUI.physicalDose)));
        assertElementsAlmostEqual(resultGUI.physicalDose,testData.resultGUI.physicalDose,'relative',1e-2,1e-2);
    
    function test_calcDoseHongPBVHEE
        testData = load('VHEE_testData.mat');
        assertTrue(DoseEngines.matRad_ParticleHongPencilBeamEngine.isAvailable(testData.pln));
        testData.pln.propDoseCalc.engine = 'HongPB';
        testData.pln.propDoseCalc.dosimetricLateralCutOff = 0.995;
        testData.pln.propDoseCalc.geometricLateralCutOff = 50;
        resultGUI = matRad_calcDoseForward(testData.ct, testData.cst, testData.stf, testData.pln, ones(sum([testData.stf(:).totalNumOfBixels]),1));
    
        assertTrue(isequal(fieldnames(resultGUI),fieldnames(testData.resultGUI)));
        assertTrue(isequal(testData.ct.cubeDim, size(resultGUI.physicalDose)));
        assertElementsAlmostEqual(resultGUI.physicalDose,testData.resultGUI.physicalDose,'relative',1e-2,1e-2);
    
    function test_calcDoseHongPBVHEE_Focused
        testData = load('VHEE_testData_Focused.mat');
        assertTrue(DoseEngines.matRad_ParticleHongPencilBeamEngine.isAvailable(testData.pln));
        testData.pln.propDoseCalc.engine = 'HongPB';
        testData.pln.propDoseCalc.dosimetricLateralCutOff = 0.995;
        testData.pln.propDoseCalc.geometricLateralCutOff = 50;
        resultGUI = matRad_calcDoseForward(testData.ct, testData.cst, testData.stf, testData.pln, ones(sum([testData.stf(:).totalNumOfBixels]),1));
    
        assertTrue(isequal(fieldnames(resultGUI),fieldnames(testData.resultGUI)));
        assertTrue(isequal(testData.ct.cubeDim, size(resultGUI.physicalDose)));
        assertElementsAlmostEqual(resultGUI.physicalDose,testData.resultGUI.physicalDose,'relative',1e-2,1e-2);
    
    function test_nonSupportedSettings
        % Radiation mode other than protons not implemented 
        testData = load('photons_testData.mat');
        testData.pln.propDoseCalc.engine = 'HongPB';    
        assertFalse(DoseEngines.matRad_ParticleHongPencilBeamEngine.isAvailable(testData.pln));
    
        % Invalid machine without radiation mode field
        testData.pln.machine = 'Empty';
        testData.pln.propDoseCalc.engine = 'HongPB';    
        assertExceptionThrown(@() DoseEngines.matRad_ParticleHongPencilBeamEngine.isAvailable(testData.pln));
        assertFalse(DoseEngines.matRad_ParticleHongPencilBeamEngine.isAvailable(testData.pln,[]));
    
    
    
    