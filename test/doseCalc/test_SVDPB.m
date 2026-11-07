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

        % The reference is computed without dij sampling (stored in testData.pln)
        testData.pln.propDoseCalc.engine = 'SVDPB';
        resultGUI = matRad_calcDoseForward(testData.ct, testData.cst, testData.stf, testData.pln, ones(sum([testData.stf(:).totalNumOfBixels]),1));

        assertTrue(isequal(fieldnames(resultGUI),fieldnames(testData.resultGUI)));
        assertTrue(isequal(testData.ct.cubeDim, size(resultGUI.physicalDose)));
        assertElementsAlmostEqual(resultGUI.physicalDose,testData.resultGUI.physicalDose,'relative',1e-2,1e-2);

    function test_dijSamplingSVDPBphotons
        testData = load('photons_testData.mat');
        testData.pln.propDoseCalc.engine = 'SVDPB';

        testData.pln.propDoseCalc.enableDijSampling = false;
        dijFull = matRad_calcDoseInfluence(testData.ct, testData.cst, testData.stf, testData.pln);
        testData.pln.propDoseCalc.enableDijSampling = true;
        dijSampled = matRad_calcDoseInfluence(testData.ct, testData.cst, testData.stf, testData.pln);
        dijSampledAgain = matRad_calcDoseInfluence(testData.ct, testData.cst, testData.stf, testData.pln);

        doseFull = dijFull.physicalDose{1};
        doseSampled = dijSampled.physicalDose{1};

        % The random number generator is seeded for each dose calculation
        assertEqual(doseSampled, dijSampledAgain.physicalDose{1});

        % Sampling only drops tail voxels, and raises the kept ones to the
        % maximum tail dose of their radiological depth cluster
        assertTrue(nnz(doseSampled) < nnz(doseFull));
        assertEqual(nnz(doseSampled & ~doseFull), 0);
        [i,j] = find(doseSampled);
        ix = sub2ind(size(doseFull),i,j);
        assertTrue(all(doseSampled(ix) >= doseFull(ix)));

        % Unbiased in expectation: the integral dose of each bixel is kept
        % within the sampling noise
        assertElementsAlmostEqual(full(sum(doseSampled,1)),full(sum(doseFull,1)),'relative',2e-2);

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
