function test_suite = test_stfGeneratorPhotonBixel

test_functions=localfunctions();

initTestSuite;

function test_basic_construct()
    stfGen = matRad_StfGeneratorPhotonSingleBeamlet();    
    assertTrue(isa(stfGen, 'matRad_StfGeneratorPhotonSingleBeamlet'));

function test_pln_construct()
    load photons_testData.mat
    stfGen = matRad_StfGeneratorPhotonSingleBeamlet(pln);
    stfGen.isAvailable(pln);
    assertTrue(isa(stfGen, 'matRad_StfGeneratorPhotonSingleBeamlet'));
    assertEqual(stfGen.gantryAngles, pln.propStf.gantryAngles);
    assertEqual(stfGen.couchAngles, pln.propStf.couchAngles);
    assertEqual(stfGen.isoCenter, pln.propStf.isoCenter);
    assertEqual(stfGen.radiationMode, pln.radiationMode);
    assertEqual(stfGen.machine, pln.machine);
    assertEqual(stfGen.bixelWidth, pln.propStf.bixelWidth);
    
function test_generate_multibeams()
    % geometry settings
    load photons_testData.mat ct cst pln;
    
    stfGen = matRad_StfGeneratorPhotonSingleBeamlet(pln);
    stf = stfGen.generate(ct,cst);
   
    assertTrue(isfield(stf, 'radiationMode'));
    assertTrue(isfield(stf, 'machine'));
    assertTrue(isfield(stf, 'gantryAngle'));
    assertTrue(isfield(stf, 'couchAngle'));
    assertTrue(isfield(stf, 'isoCenter'));
    assertTrue(isfield(stf, 'bixelWidth'));
    assertTrue(isfield(stf, 'SAD'));
    assertTrue(isfield(stf, 'SCD'));
    assertTrue(isfield(stf, 'numOfRays'));
    assertTrue(isfield(stf, 'numOfBixelsPerRay'));
    assertTrue(isfield(stf, 'totalNumOfBixels'));
    assertTrue(isfield(stf, 'sourcePoint'));
    assertTrue(isfield(stf, 'sourcePoint_bev'));
    assertTrue(isfield(stf, 'ray'));

    for i = 1:numel(stf)
        
        assertEqual(stf(i).totalNumOfBixels,1);
        assertEqual(stf(i).numOfBixelsPerRay,1);
        assertEqual(stf(i).numOfRays,1);
        assertEqual(stf(i).bixelWidth,stfGen.bixelWidth);
        assertEqual(stf(i).radiationMode,stfGen.radiationMode);
        assertEqual(stf(i).machine,pln.machine);
        assertEqual(stf(i).gantryAngle,stfGen.gantryAngles(i));
        assertEqual(stf(i).couchAngle,stfGen.couchAngles(i));

        rotMat = matRad_getRotationMatrix(stf(i).gantryAngle,stf(i).couchAngle);
        assertEqual(stf(i).sourcePoint,stf(i).sourcePoint_bev*rotMat);

        assertEqual(stf(i).sourcePoint_bev,[0 -stf(i).SAD 0]);        

        assertTrue(isstruct(stf(i).ray) && numel(stf(i).ray) == 1);
        assertEqual(stf(i).ray.rayPos, [0 0 0]);
        assertEqual(stf(i).ray.rayPos_bev, [0 0 0]);
        assertEqual(stf(i).ray.targetPoint_bev, [0 stf(i).SAD 0]);
        assertEqual(stf(i).ray.targetPoint, stf(i).ray.targetPoint_bev*rotMat);

        assertTrue(isfield(stf(i).ray,'beamletCornersAtIso'));
        assertTrue(isfield(stf(i).ray,'rayCorners_SCD'));
        assertTrue(isscalar(stf(i).ray.energy));
    end
        
    function test_generate_single_beam()
        % geometry settings
        load photons_testData.mat ct cst pln;
        
        stfGen = matRad_StfGeneratorPhotonSingleBeamlet(pln);

        stfGen.gantryAngles = 0;
        % couch angles are not auto-synced; the mismatch is caught on use
        assertExceptionThrown(@() stfGen.numOfBeams, 'matRad:Error');
        stfGen.couchAngles = 0;
        assertEqual(stfGen.numOfBeams, 1);

        stf = stfGen.generate(ct,cst);
    
        assertTrue(isfield(stf, 'radiationMode'));
        assertTrue(isfield(stf, 'machine'));
        assertTrue(isfield(stf, 'gantryAngle'));
        assertTrue(isfield(stf, 'couchAngle'));
        assertTrue(isfield(stf, 'isoCenter'));
        assertTrue(isfield(stf, 'bixelWidth'));
        assertTrue(isfield(stf, 'SAD'));
        assertTrue(isfield(stf, 'SCD'));
        assertTrue(isfield(stf, 'numOfRays'));
        assertTrue(isfield(stf, 'numOfBixelsPerRay'));
        assertTrue(isfield(stf, 'totalNumOfBixels'));
        assertTrue(isfield(stf, 'sourcePoint'));
        assertTrue(isfield(stf, 'sourcePoint_bev'));
        assertTrue(isfield(stf, 'ray'));
    
        assertEqual(stf.totalNumOfBixels,1);
        assertEqual(stf.numOfBixelsPerRay,1);
        assertEqual(stf.numOfRays,1);
        assertEqual(stf.bixelWidth,stfGen.bixelWidth);
        assertEqual(stf.radiationMode,stfGen.radiationMode);
        assertEqual(stf.machine,pln.machine);
        assertEqual(stf.gantryAngle,stfGen.gantryAngles);
        assertEqual(stf.couchAngle,stfGen.couchAngles);

        rotMat = matRad_getRotationMatrix(stf.gantryAngle,stf.couchAngle);
        assertEqual(stf.sourcePoint_bev,[0 -stf.SAD 0]);        
        assertEqual(stf.sourcePoint,stf.sourcePoint_bev*rotMat);

        assertTrue(isstruct(stf.ray) && numel(stf.ray) == 1);
        assertEqual(stf.ray.rayPos, [0 0 0]);
        assertEqual(stf.ray.rayPos_bev, [0 0 0]);
        assertEqual(stf.ray.targetPoint_bev, [0 stf.SAD 0]);
        assertEqual(stf.ray.targetPoint, stf.ray.targetPoint_bev*rotMat);
        assertTrue(isfield(stf.ray,'beamletCornersAtIso'));
        assertTrue(isfield(stf.ray,'rayCorners_SCD'));

function test_angle_consistency_pln()
    load photons_testData.mat pln;

    % Inconsistent numbers of gantry and couch angles in pln throw an error
    plnBad = pln;
    plnBad.propStf.gantryAngles = [0 90 270];
    plnBad.propStf.couchAngles  = [0 0];
    assertExceptionThrown(@() matRad_StfGeneratorPhotonSingleBeamlet(plnBad), 'matRad:Error');

    % A scalar couch angle is valid for multiple gantry angles
    plnScalar = pln;
    plnScalar.propStf.gantryAngles = [0 90 270];
    plnScalar.propStf.couchAngles  = 0;
    stfGen = matRad_StfGeneratorPhotonSingleBeamlet(plnScalar);
    assertEqual(stfGen.numOfBeams, 3);

function test_angle_consistency_properties()
    stfGen = matRad_StfGeneratorPhotonSingleBeamlet();

    % Sequential assignment with transient mismatch is allowed
    stfGen.gantryAngles = [0 90 270];
    stfGen.couchAngles  = [0 10 20];
    assertEqual(stfGen.numOfBeams, 3);

    % Mismatch is caught when the number of beams is queried
    stfGen.couchAngles = [0 10];
    assertExceptionThrown(@() stfGen.numOfBeams, 'matRad:Error');

function test_scalar_couch_angle_expansion()
    load photons_testData.mat ct cst pln;

    stfGen = matRad_StfGeneratorPhotonSingleBeamlet(pln);
    stfGen.gantryAngles = [0 180];
    stfGen.couchAngles  = 10;

    stf = stfGen.generate(ct,cst);

    assertEqual([stf.gantryAngle], [0 180]);
    assertEqual([stf.couchAngle], [10 10]);
    assertEqual(stfGen.couchAngles, [10 10]);