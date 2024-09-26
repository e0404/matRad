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
            assertTrue(numel(stfGen.couchAngles) == 1);
            stfGen.couchAngles = 0;
    
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