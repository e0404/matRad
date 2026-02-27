function test_suite = test_FSPB

    test_functions=localfunctions();
    
    initTestSuite;
    
    function test_getSubsamplingPBEngineFromPln
        % Single gaussian lateral model
        testData.pln = struct('radiationMode','protons','machine','Generic');
        testData.pln.propDoseCalc.engine = 'SubsamplingPB';
        engine = DoseEngines.matRad_ParticleFineSamplingPencilBeamEngine.getEngineFromPln(testData.pln);
        assertTrue(isa(engine,'DoseEngines.matRad_ParticleFineSamplingPencilBeamEngine'));
        
    function test_loadMachineForSubsamplingPB
        possibleRadModes = DoseEngines.matRad_ParticleFineSamplingPencilBeamEngine.possibleRadiationModes;
        for i = 1:numel(possibleRadModes)
            machine = DoseEngines.matRad_ParticleFineSamplingPencilBeamEngine.loadMachine(possibleRadModes{i},'Generic');
            assertTrue(isstruct(machine));
            assertTrue(isfield(machine, 'meta'));
            assertTrue(isfield(machine.meta, 'radiationMode'));
            assertTrue(strcmp(machine.meta.radiationMode, possibleRadModes{i}));
        end

    
    function test_calcDoseSubsamplingPBprotonsFitCircle
        testData = load('protons_testData.mat');

        assertTrue(DoseEngines.matRad_ParticleFineSamplingPencilBeamEngine.isAvailable(testData.pln));
        
        testData.pln.propDoseCalc.engine = 'SubsamplingPB';
        testData.pln.propDoseCalc.dosimetricLateralCutOff = 0.995;
        testData.pln.propDoseCalc.geometricLateralCutOff = 50;
        testData.pln.propDoseCalc.fineSampling.method = 'fitCircle';
        
        for N = [2,3,8]
            testData.pln.propDoseCalc.fineSampling.N = N;
            
            resultGUI = matRad_calcDoseForward(testData.ct, testData.cst, testData.stf, testData.pln, ones(sum([testData.stf(:).totalNumOfBixels]),1));
    
            assertTrue(isequal(fieldnames(resultGUI),fieldnames(testData.resultGUI)));
            assertTrue(isequal(testData.ct.cubeDim, size(resultGUI.physicalDose)));
            assertElementsAlmostEqual(resultGUI.physicalDose,testData.resultGUI.physicalDose,'relative',1e-2,1e-2);
        end
    
    function test_calcDoseSubsamplingPBprotonsFitSquare
        testData = load('protons_testData.mat');

        assertTrue(DoseEngines.matRad_ParticleFineSamplingPencilBeamEngine.isAvailable(testData.pln));
        
        testData.pln.propDoseCalc.engine = 'SubsamplingPB';
        testData.pln.propDoseCalc.dosimetricLateralCutOff = 0.995;
        testData.pln.propDoseCalc.geometricLateralCutOff = 50;
        testData.pln.propDoseCalc.fineSampling.method = 'fitSquare';
        
        for N = [2,3]
            testData.pln.propDoseCalc.fineSampling.N = N;
            resultGUI = matRad_calcDoseForward(testData.ct, testData.cst, testData.stf, testData.pln, ones(sum([testData.stf(:).totalNumOfBixels]),1));
    
            assertTrue(isequal(fieldnames(resultGUI),fieldnames(testData.resultGUI)));
            assertTrue(isequal(testData.ct.cubeDim, size(resultGUI.physicalDose)));
            assertElementsAlmostEqual(resultGUI.physicalDose,testData.resultGUI.physicalDose,'relative',1e-2,1e-2);
        end

    function test_calcDoseSubsamplingPBprotonsRusso
        testData = load('protons_testData.mat');
        
        testData.pln.propDoseCalc.fineSampling.N = 10;
        testData.pln.propDoseCalc.fineSampling.sigmaSub = 2;
        testData.pln.propDoseCalc.fineSampling.method = 'russo';

        resultGUI = matRad_calcDoseForward(testData.ct, testData.cst, testData.stf, testData.pln, ones(sum([testData.stf(:).totalNumOfBixels]),1));
    
        assertTrue(isequal(fieldnames(resultGUI),fieldnames(testData.resultGUI)));
        assertTrue(isequal(testData.ct.cubeDim, size(resultGUI.physicalDose)));
        assertElementsAlmostEqual(resultGUI.physicalDose,testData.resultGUI.physicalDose,'relative',1e-2,1e-2);

    function test_calcDoseSubsamplingPBhelium
        testData = load('helium_testData.mat');
        assertTrue(DoseEngines.matRad_ParticleFineSamplingPencilBeamEngine.isAvailable(testData.pln));

        testData.pln.propDoseCalc.engine = 'SubsamplingPB';
        testData.pln.propDoseCalc.dosimetricLateralCutOff = 0.995;
        testData.pln.propDoseCalc.geometricLateralCutOff = 50;
        testData.pln.propDoseCalc.fineSampling.N = 10;
        testData.pln.propDoseCalc.fineSampling.sigmaSub = 2;
        testData.pln.propDoseCalc.fineSampling.method = 'russo';
        resultGUI = matRad_calcDoseForward(testData.ct, testData.cst, testData.stf, testData.pln, ones(sum([testData.stf(:).totalNumOfBixels]),1));
    
        assertTrue(isequal(fieldnames(resultGUI),fieldnames(testData.resultGUI)));
        assertTrue(isequal(testData.ct.cubeDim, size(resultGUI.physicalDose)));
        assertElementsAlmostEqual(resultGUI.physicalDose,testData.resultGUI.physicalDose,'relative',1e-2,1e-2);

    function test_calcDoseSubsamplingPBcarbon
        testData = load('carbon_testData.mat');
        assertTrue(DoseEngines.matRad_ParticleFineSamplingPencilBeamEngine.isAvailable(testData.pln));

        testData.pln.propDoseCalc.engine = 'SubsamplingPB';
        testData.pln.propDoseCalc.dosimetricLateralCutOff = 0.995;
        testData.pln.propDoseCalc.geometricLateralCutOff = 50;
        testData.pln.propDoseCalc.fineSampling.N = 10;
        testData.pln.propDoseCalc.fineSampling.sigmaSub = 2;
        testData.pln.propDoseCalc.fineSampling.method = 'russo';        
        resultGUI = matRad_calcDoseForward(testData.ct, testData.cst, testData.stf, testData.pln, ones(sum([testData.stf(:).totalNumOfBixels]),1));
    
        assertTrue(isequal(fieldnames(resultGUI),fieldnames(testData.resultGUI)));
        assertTrue(isequal(testData.ct.cubeDim, size(resultGUI.physicalDose)));
        assertElementsAlmostEqual(resultGUI.physicalDose,testData.resultGUI.physicalDose,'relative',1e-2,1e-2);
   
    
    function test_nonSupportedSettings
        % Radiation mode other than protons not implemented 
        testData = load('photons_testData.mat');
        testData.pln.propDoseCalc.engine = 'SubsamplingPB';    
        assertFalse(DoseEngines.matRad_ParticleFineSamplingPencilBeamEngine.isAvailable(testData.pln));
    
        % Invalid machine without radiation mode field
        testData.pln.machine = 'Empty';
        testData.pln.propDoseCalc.engine = 'SubsamplingPB';    
        assertExceptionThrown(@() DoseEngines.matRad_ParticleFineSamplingPencilBeamEngine.isAvailable(testData.pln));
        assertFalse(DoseEngines.matRad_ParticleFineSamplingPencilBeamEngine.isAvailable(testData.pln,[]));
    
    
    
    