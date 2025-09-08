function test_suite = test_stfGeneratorVHEE

    test_functions=localfunctions();
    
    initTestSuite;
    
    function test_basic_construct()
        stfGen = matRad_StfGeneratorParticleVHEE();    
        assertTrue(isa(stfGen, 'matRad_StfGeneratorParticleVHEE'));
    
    function test_pln_construct()
        load VHEE_testData.mat
        pln.propStf.energy = 150;
        stfGen = matRad_StfGeneratorParticleVHEE(pln);
        stfGen.isAvailable(pln);
        assertTrue(isa(stfGen, 'matRad_StfGeneratorParticleVHEE'));
        assertEqual(stfGen.gantryAngles, pln.propStf.gantryAngles);
        assertEqual(stfGen.couchAngles, pln.propStf.couchAngles);
        assertEqual(stfGen.isoCenter, pln.propStf.isoCenter);
        assertEqual(stfGen.radiationMode, pln.radiationMode);
        assertEqual(stfGen.machine, pln.machine);
        assertEqual(stfGen.bixelWidth, pln.propStf.bixelWidth);
        assertEqual(stfGen.energy, pln.propStf.energy);

    function test_pln_construct_focused()
        load VHEE_testData.mat
        pln.machine = 'Focused';
        pln.propStf.energy = 150;
        stfGen = matRad_StfGeneratorParticleVHEE(pln);
        stfGen.isAvailable(pln);
        assertTrue(isa(stfGen, 'matRad_StfGeneratorParticleVHEE'));
        assertEqual(stfGen.gantryAngles, pln.propStf.gantryAngles);
        assertEqual(stfGen.couchAngles, pln.propStf.couchAngles);
        assertEqual(stfGen.isoCenter, pln.propStf.isoCenter);
        assertEqual(stfGen.radiationMode, pln.radiationMode);
        assertEqual(stfGen.machine, pln.machine);
        assertEqual(stfGen.bixelWidth, pln.propStf.bixelWidth);
        assertEqual(stfGen.energy, pln.propStf.energy);

    function test_generate_multibeams()
        % geometry settings
        load VHEE_testData.mat ct cst pln stf;
        
        stfGen = matRad_StfGeneratorParticleVHEE(pln);
        stf2 = stfGen.generate(ct,cst);
       
        assertTrue(isfield(stf2, 'radiationMode'));
        assertTrue(isfield(stf2, 'machine'));
        assertTrue(isfield(stf2, 'gantryAngle'));
        assertTrue(isfield(stf2, 'couchAngle'));
        assertTrue(isfield(stf2, 'isoCenter'));
        assertTrue(isfield(stf2, 'bixelWidth'));
        assertTrue(isfield(stf2, 'SAD'));
        assertTrue(isfield(stf2, 'numOfRays'));
        assertTrue(isfield(stf2, 'numOfBixelsPerRay'));
        assertTrue(isfield(stf2, 'totalNumOfBixels'));
        assertTrue(isfield(stf2, 'sourcePoint'));
        assertTrue(isfield(stf2, 'sourcePoint_bev'));
        assertTrue(isfield(stf2, 'ray'));
    
        for i = 1:numel(stf2)
            
            assertEqual(stf2(i).totalNumOfBixels,stf(i).totalNumOfBixels);
            assertEqual(stf2(i).numOfBixelsPerRay,stf(i).numOfBixelsPerRay);
            assertEqual(stf2(i).numOfRays,stf(i).numOfRays);
            assertEqual(stf2(i).bixelWidth,stfGen.bixelWidth);
            assertEqual(stf2(i).radiationMode,stfGen.radiationMode);
            assertEqual(stf2(i).machine,pln.machine);
            assertEqual(stf2(i).gantryAngle,stfGen.gantryAngles(i));
            assertEqual(stf2(i).couchAngle,stfGen.couchAngles(i));
            assertTrue(all([stf2(i).ray.energy] == stfGen.energy));
    
            rotMat = matRad_getRotationMatrix(stf2(i).gantryAngle,stf2(i).couchAngle);
            assertEqual(stf2(i).sourcePoint,stf2(i).sourcePoint_bev*rotMat);
            assertEqual(stf2(i).sourcePoint_bev,[0 -stf2(i).SAD 0]);        
    
            assertTrue(isstruct(stf2(i).ray));
            assertEqual(numel(stf2(i).ray),numel(stf(i).ray));
            assertEqual(numel(stf2(i).ray),stf2(i).numOfRays);

            rayPosTest = vertcat(stf2(i).ray.rayPos);
            rayPosTest_bev = rayPosTest*rotMat;
            rayPos_bevTest = vertcat(stf2(i).ray.rayPos_bev);
            rayPosRef = vertcat(stf(i).ray.rayPos);
            assertElementsAlmostEqual(sort(rayPosTest,1),sort(rayPosRef,1));
            assertElementsAlmostEqual(sort(rayPos_bevTest,1),sort(rayPosTest_bev,1));

            targetPointTest = vertcat(stf2(i).ray.targetPoint);
            targetPointRef = vertcat(stf(i).ray.targetPoint);
            targetPoint_bevTest = vertcat(stf2(i).ray.targetPoint_bev);
            targetPoint_bevRef = vertcat(stf(i).ray.targetPoint_bev);
            assertElementsAlmostEqual(sort(targetPointTest,1),sort(targetPointRef,1));
            assertElementsAlmostEqual(sort(targetPoint_bevTest,1),sort(targetPoint_bevRef,1));
            assertElementsAlmostEqual(sort(targetPoint_bevTest,1),sort(targetPointTest*rotMat,1));
            
            energiesTest = [stf2(i).ray.energy];
            energiesRef = [stf(i).ray.energy];
            assertEqual(unique(energiesTest),unique(energiesRef));
            
        end
