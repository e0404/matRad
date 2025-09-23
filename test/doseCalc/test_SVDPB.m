function test_suite = test_SVDPB

    test_functions=localfunctions();

    initTestSuite;

    function test_getSVDPBEngineFromPln
        % Single gaussian lateral model
        testData.pln = struct('radiationMode','photons','machine','Generic');
        testData.pln.propDoseCalc.engine = 'SVDPB';
        engine = DoseEngines.matRad_PhotonPencilBeamSVDEngine.getEngineFromPln(testData.pln);
        assertTrue(isa(engine,'DoseEngines.matRad_PhotonPencilBeamSVDEngine'));

    function test_loadMachineForSVDPB
        possibleRadModes = DoseEngines.matRad_PhotonPencilBeamSVDEngine.possibleRadiationModes;
        for i = 1:numel(possibleRadModes)
            machine = DoseEngines.matRad_PhotonPencilBeamSVDEngine.loadMachine(possibleRadModes{i},'Generic');
            assertTrue(isstruct(machine));
            assertTrue(isfield(machine, 'meta'));
            assertTrue(isfield(machine.meta, 'radiationMode'));
            assertTrue(strcmp(machine.meta.radiationMode, possibleRadModes{i}));
        end

    function test_calcDoseSVDPBphotons
        testData = load('photons_testData.mat');

        assertTrue(DoseEngines.matRad_PhotonPencilBeamSVDEngine.isAvailable(testData.pln));

        testData.pln.propDoseCalc.engine = 'SVDPB';
        testData.pln.propDoseCalc.dosimetricLateralCutOff = 0.995;
        testData.pln.propDoseCalc.geometricLateralCutOff = 50;
        testData.pln.propDoseCalc.kernelCutOff = Inf;
        if moxunit_util_platform_is_octave()
          %The random number generator is not consistent between octave and matlab
          testData.pln.propDoseCalc.enableDijSampling = false;
        end
        resultGUI = matRad_calcDoseForward(testData.ct, testData.cst, testData.stf, testData.pln, ones(sum([testData.stf(:).totalNumOfBixels]),1));

        assertTrue(isequal(fieldnames(resultGUI),fieldnames(testData.resultGUI)));
        assertTrue(isequal(testData.ct.cubeDim, size(resultGUI.physicalDose)));
        assertElementsAlmostEqual(resultGUI.physicalDose,testData.resultGUI.physicalDose,'relative',1e-2,1e-2);


    function test_nonSupportedSettings
        % Radiation mode other than photons not implemented
        testData = load('protons_testData.mat');
        testData.pln.propDoseCalc.engine = 'SVDPB';
        assertFalse(DoseEngines.matRad_PhotonPencilBeamSVDEngine.isAvailable(testData.pln));

        % Invalid machine without radiation mode field
        testData.pln.machine = 'Empty';
        testData.pln.propDoseCalc.engine = 'SVDPB';
        assertExceptionThrown(@() DoseEngines.matRad_PhotonPencilBeamSVDEngine.isAvailable(testData.pln));
        assertFalse(DoseEngines.matRad_PhotonPencilBeamSVDEngine.isAvailable(testData.pln,[]));




