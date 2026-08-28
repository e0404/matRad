function test_suite = test_stfGeneratorParticleBeamlet

    test_functions=localfunctions();
    
    initTestSuite;
    
    function test_basic_construct()
        stfGen = matRad_StfGeneratorParticleSingleBeamlet();    
        assertTrue(isa(stfGen, 'matRad_StfGeneratorParticleSingleBeamlet'));
        assertEqual(stfGen.radiationMode,'protons');
    
    function test_pln_construct()
        load protons_testData.mat
        stfGen = matRad_StfGeneratorParticleSingleBeamlet(pln);
        assertTrue(stfGen.isAvailable(pln));
        assertTrue(isa(stfGen, 'matRad_StfGeneratorParticleSingleBeamlet'));
        assertEqual(stfGen.gantryAngles, pln.propStf.gantryAngles);
        assertEqual(stfGen.couchAngles, pln.propStf.couchAngles);
        assertEqual(stfGen.isoCenter, pln.propStf.isoCenter);
        assertEqual(stfGen.radiationMode, pln.radiationMode);
        assertEqual(stfGen.machine, pln.machine);
        assertEqual(stfGen.bixelWidth, pln.propStf.bixelWidth);
        assertEqual(stfGen.radiationMode,'protons');

        pln.radiationMode = 'helium';
        stfGen = matRad_StfGeneratorParticleSingleBeamlet(pln);
        assertTrue(stfGen.isAvailable(pln));
        assertTrue(isa(stfGen, 'matRad_StfGeneratorParticleSingleBeamlet'));
        assertEqual(stfGen.radiationMode, pln.radiationMode);

        pln.radiationMode = 'carbon';
        stfGen = matRad_StfGeneratorParticleSingleBeamlet(pln);
        assertTrue(stfGen.isAvailable(pln));
        assertTrue(isa(stfGen, 'matRad_StfGeneratorParticleSingleBeamlet'));
        assertEqual(stfGen.radiationMode, pln.radiationMode);        

        
    function test_generate_multibeams()
        % geometry settings
        load protons_testData.mat ct cst pln;
        
        stfGen = matRad_StfGeneratorParticleSingleBeamlet(pln);
        stf = stfGen.generate(ct,cst);
       
    
        assertTrue(isfield(stf, 'radiationMode'));
        assertTrue(isfield(stf, 'machine'));
        assertTrue(isfield(stf, 'gantryAngle'));
        assertTrue(isfield(stf, 'couchAngle'));
        assertTrue(isfield(stf, 'isoCenter'));
        assertTrue(isfield(stf, 'bixelWidth'));
        assertTrue(isfield(stf, 'SAD'));
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
    
            assertTrue(isfield(stf(i).ray,'rangeShifter'));            
            assertTrue(isscalar(stf(i).ray.energy));
        end
            
        function test_generate_single_beam()
            % geometry settings
            load protons_testData.mat ct cst pln;
            
            stfGen = matRad_StfGeneratorParticleSingleBeamlet(pln);
    
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
            
            assertTrue(isfield(stf.ray,'rangeShifter'));
            assertTrue(isscalar(stf.ray.energy));

function test_construct_oxygen()
    load protons_testData.mat pln;
    pln.radiationMode = 'oxygen';
    stfGen = matRad_StfGeneratorParticleSingleBeamlet(pln);
    assertTrue(stfGen.isAvailable(pln));
    assertTrue(isa(stfGen, 'matRad_StfGeneratorParticleSingleBeamlet'));
    assertEqual(stfGen.radiationMode, 'oxygen');

function test_generate_explicitEnergy()
    % specifying the energy directly picks the closest available energy
    load protons_testData.mat ct cst pln;
    stfGen = matRad_StfGeneratorParticleSingleBeamlet(pln);
    stfGen.energy = 100;
    stf = stfGen.generate(ct,cst);
    assertTrue(isscalar(stf(1).ray.energy));
    % closest available energy to 100 MeV
    assertElementsAlmostEqual(stf(1).ray.energy, 99.7909, 'relative', 1e-3);

function test_generate_autoIsoCenter()
    % an empty isoCenter must be computed automatically from the cst
    load protons_testData.mat ct cst pln;
    stfGen = matRad_StfGeneratorParticleSingleBeamlet(pln);
    stfGen.isoCenter = [];
    stf = stfGen.generate(ct,cst);
    assertTrue(isstruct(stf));
    assertTrue(isfield(stf, 'isoCenter'));

function test_generate_rangeShifter()
    % range shifter placement fills the rangeShifter sub-struct
    load protons_testData.mat ct cst pln;
    stfGen = matRad_StfGeneratorParticleSingleBeamlet(pln);
    stfGen.useRangeShifter = true;
    stf = stfGen.generate(ct,cst);
    assertEqual(stf(1).ray.rangeShifter.ID, 1);
    assertEqual(stf(1).ray.rangeShifter.eqThickness, stfGen.raShiThickness);

function test_isAvailable_invalidMachine()
    % an invalid machine struct (no meta/data) must be reported unavailable
    load protons_testData.mat pln;
    [available, msg] = matRad_StfGeneratorParticleSingleBeamlet.isAvailable(pln, struct('foo',1));
    assertFalse(available);
    assertTrue(ischar(msg));