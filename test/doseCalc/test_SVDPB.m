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
        % Summed over all bixels, the sampling noise is far below 0.15 %
        assertElementsAlmostEqual(full(sum(doseSampled(:))),full(sum(doseFull(:))),'relative',1.5e-3);

    function test_dijSamplingShallowRadDepths
        % All tail voxels lie shallower than two radiological depth
        % clusters; sampling must still keep part of the tail
        ct = struct();
        ct.cubeDim = [3 51 51];
        ct.resolution = struct('x',2,'y',2,'z',2);
        ct.cubeHU = {zeros(ct.cubeDim)};
        ct.numOfCtScen = 1;
        cst = {0,'Body','OAR',{(1:prod(ct.cubeDim))'},struct('Priority',2),[]; ...
               1,'Target','TARGET',{sub2ind(ct.cubeDim,2,26,26)},struct('Priority',1),[]};

        pln = struct('radiationMode','photons','machine','Generic','numOfFractions',1);
        pln.propStf.gantryAngles = 0;
        pln.propStf.couchAngles = 0;
        pln.propStf.bixelWidth = 5;
        pln.propStf.isoCenter = matRad_getIsoCenter(cst,ct,0);
        pln.propDoseCalc.engine = 'SVDPB';
        pln.propDoseCalc.doseGrid.resolution = ct.resolution;
        pln.propDoseCalc.geometricLateralCutOff = 50;
        pln.propDoseCalc.kernelCutOff = 50;
        pln.multScen = 'nomScen';

        stf = matRad_generateStf(ct,cst,pln);
        pln.propDoseCalc.enableDijSampling = false;
        dijFull = matRad_calcDoseInfluence(ct,cst,stf,pln);
        pln.propDoseCalc.enableDijSampling = true;
        dijSampled = matRad_calcDoseInfluence(ct,cst,stf,pln);

        assertElementsAlmostEqual(full(sum(dijSampled.physicalDose{1}(:))),full(sum(dijFull.physicalDose{1}(:))),'relative',1e-2);

    function test_bixelDoseMirrorSymmetricWithKernelCutOff
        % The Gauss filter and kernel grids must be symmetric around the
        % ray, also when the kernel is cut off where it is not yet zero
        [ct, cst, pln, stf] = helper_centralRaySetup();
        pln.propDoseCalc.kernelCutOff = 20;

        resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, 1);
        dose = resultGUI.physicalDose;
        assertElementsAlmostEqual(dose, flip(dose, 2), 'absolute', 1e-5 * max(dose(:)));
        assertElementsAlmostEqual(dose, flip(dose, 3), 'absolute', 1e-5 * max(dose(:)));

    function test_fieldBasedDoseMirrorSymmetric
        % A field shape pixel covers [g, g + convResolution), so a symmetric
        % open field must give a mirror symmetric dose
        [ct, cst, pln, stf] = helper_centralRaySetup();
        [pln, stf] = helper_toFieldBased(pln, stf, ones(20));

        resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, 1);
        dose = resultGUI.physicalDose;
        assertElementsAlmostEqual(dose, flip(dose, 2), 'absolute', 1e-5 * max(dose(:)));
        assertElementsAlmostEqual(dose, flip(dose, 3), 'absolute', 1e-5 * max(dose(:)));

    function test_fieldBasedShapeSizeMismatchErrors
        % fieldWidth = 10 mm at 0.5 mm resolution requires 20 x 20 pixels
        [ct, cst, pln, stf] = helper_centralRaySetup();
        [pln, stf] = helper_toFieldBased(pln, stf, ones(21));

        assertExceptionThrown(@() matRad_calcDoseForward(ct, cst, stf, pln, 1), 'matRad:Error');

    function test_fieldBasedDoseMatchesSequencedBixelDose
        % Field-based dose of the sequenced apertures must reproduce the
        % bixel-based dose of the sequenced weights
        p = load('photons_testData.mat');
        pln = p.pln;
        pln.propDoseCalc.engine = 'SVDPB';
        pln.propSeq.sequencer = 'siochi';
        resultGUI = matRad_sequencing(p.resultGUI, p.stf, pln);
        sequencing = resultGUI.sequencing;

        resultBixel = matRad_calcDoseForward(p.ct, p.cst, p.stf, pln, sequencing.w);
        [plnField, stfField] = matRad_aperture2collimation(pln, p.stf, sequencing, sequencing.apertureInfo);
        resultField = matRad_calcDoseForward(p.ct, p.cst, stfField, plnField);

        doseBixel = resultBixel.physicalDose;
        assertElementsAlmostEqual(resultField.physicalDose, doseBixel, 'absolute', 5e-3 * max(doseBixel(:)));

    function [pln, stf] = helper_toFieldBased(pln, stf, shape)
        fieldWidth = 10;
        pln.propStf.bixelWidth = 'field';
        pln.propStf.collimation = struct('convResolution', 0.5, 'fieldWidth', fieldWidth);
        stf.bixelWidth = 'field';
        stf.ray.shape = shape;

    function [ct, cst, pln, stf] = helper_centralRaySetup()
        % Water cube with the isocenter in its center voxel and a single
        % 5 mm bixel on the central axis
        ct = struct();
        ct.cubeDim = [41 41 41];
        ct.resolution = struct('x', 1, 'y', 1, 'z', 1);
        ct.cubeHU = {zeros(ct.cubeDim)};
        ct.numOfCtScen = 1;
        cst = {0, 'Body', 'OAR', {(1:prod(ct.cubeDim))'}, struct('Priority', 2), []
               1, 'Target', 'TARGET', {sub2ind(ct.cubeDim, 21, 21, 21)}, struct('Priority', 1), []};

        pln = struct('radiationMode', 'photons', 'machine', 'Generic', 'numOfFractions', 1);
        pln.propStf.gantryAngles = 0;
        pln.propStf.couchAngles = 0;
        pln.propStf.bixelWidth = 5;
        pln.propStf.isoCenter = matRad_getIsoCenter(cst, ct, 0);
        pln.propDoseCalc.engine = 'SVDPB';
        pln.propDoseCalc.enableDijSampling = false;
        pln.propDoseCalc.doseGrid.resolution = ct.resolution;
        pln.propDoseCalc.geometricLateralCutOff = 20;
        pln.multScen = 'nomScen';

        stf = matRad_generateStf(ct, cst, pln);
        [~, ix] = min(sum(reshape([stf.ray.rayPos_bev], 3, []).^2, 1));
        stf.ray = stf.ray(ix);
        stf.numOfRays = 1;
        stf.numOfBixelsPerRay = 1;
        stf.totalNumOfBixels = 1;

    function test_customPrimaryFluenceGantryInvariance
        % The primary fluence of an off-axis bixel must not depend on the
        % gantry angle in a phantom that is symmetric under gantry rotation
        restoreUserFolders = helper_addRampFluenceMachine(); %#ok<NASGU>

        gantryAngles = [0 90 180 270];
        integralDose = zeros(size(gantryAngles));
        for i = 1:numel(gantryAngles)
            dose = helper_offAxisBixelDose(gantryAngles(i),[30 0 0],true);
            integralDose(i) = sum(dose(:));
        end
        assertElementsAlmostEqual(integralDose,integralDose(1)*ones(size(gantryAngles)),'relative',1e-3);

    function test_customPrimaryFluenceTiltDirection
        % The primary fluence grows with the distance to the central axis,
        % so it pulls the dose of an off-axis bixel outward
        restoreUserFolders = helper_addRampFluenceMachine(); %#ok<NASGU>

        for rayPosX = [-30 30]
            [doseCustom,x] = helper_offAxisBixelDose(0,[rayPosX 0 0],true);
            doseUniform = helper_offAxisBixelDose(0,[rayPosX 0 0],false);

            centroidCustom = sum(x(:).*doseCustom(:))/sum(doseCustom(:));
            centroidUniform = sum(x(:).*doseUniform(:))/sum(doseUniform(:));
            assertTrue(sign(centroidCustom - centroidUniform) == sign(rayPosX));
        end

    function restoreUserFolders = helper_addRampFluenceMachine()
        % Adds a temporary user folder with a photon machine whose primary
        % fluence rises linearly with the distance to the central axis
        matRad_cfg = MatRad_Config.instance();
        oldUserFolders = matRad_cfg.userfolders;
        tmpFolder = helper_temporaryFolder('matRad_test_SVDPB_rampFluence');
        mkdir(fullfile(tmpFolder,'machines'));

        machine = matRad_loadMachine(struct('radiationMode','photons','machine','Generic'));
        machine.meta.name = 'RampFluence';
        r = (0:400)';
        machine.data.primaryFluence = [r r];
        save(fullfile(tmpFolder,'machines','photons_RampFluence.mat'),'machine');

        matRad_cfg.userfolders = [oldUserFolders(:)' {tmpFolder}];
        restoreUserFolders = onCleanup(@() helper_restoreUserFolders(oldUserFolders));

    function helper_restoreUserFolders(userFolders)
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.userfolders = userFolders;

    function [dose,x] = helper_offAxisBixelDose(gantryAngle,rayPosBev,useCustomPrimaryFluence)
        % Dose of the single bixel at rayPosBev [mm] in a water cube whose
        % center is the isocenter, and the world x coordinate of each voxel
        ct = struct();
        ct.cubeDim = [51 51 51];
        ct.resolution = struct('x',2,'y',2,'z',2);
        ct.cubeHU = {zeros(ct.cubeDim)};
        ct.numOfCtScen = 1;

        [yy,xx,zz] = ndgrid(1:ct.cubeDim(1),1:ct.cubeDim(2),1:ct.cubeDim(3));
        isTarget = (xx-26).^2 + (yy-26).^2 + (zz-26).^2 <= 17^2;
        cst = {0,'Body','OAR',{(1:prod(ct.cubeDim))'},struct('Priority',2),[]; ...
               1,'Target','TARGET',{find(isTarget)},struct('Priority',1),[]};

        pln = struct('radiationMode','photons','machine','RampFluence','numOfFractions',1);
        pln.propStf.gantryAngles = gantryAngle;
        pln.propStf.couchAngles = 0;
        pln.propStf.bixelWidth = 5;
        pln.propStf.isoCenter = matRad_getIsoCenter(cst,ct,0);
        pln.propDoseCalc.engine = 'SVDPB';
        pln.propDoseCalc.enableDijSampling = false;
        pln.propDoseCalc.useCustomPrimaryPhotonFluence = useCustomPrimaryFluence;
        pln.propDoseCalc.doseGrid.resolution = ct.resolution;
        pln.multScen = 'nomScen';

        stf = matRad_generateStf(ct,cst,pln);
        rayPos = reshape([stf.ray.rayPos_bev],3,[])';
        ix = find(all(abs(rayPos - rayPosBev) < 1e-6,2));
        assertEqual(numel(ix),1);
        stf.ray = stf.ray(ix);
        stf.numOfRays = 1;
        stf.numOfBixelsPerRay = 1;
        stf.totalNumOfBixels = 1;

        resultGUI = matRad_calcDoseForward(ct,cst,stf,pln,1);
        dose = resultGUI.physicalDose;
        x = (xx - 26)*ct.resolution.x;

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
