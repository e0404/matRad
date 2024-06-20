function test_suite = test_nominalScenario

test_functions=localfunctions();

initTestSuite;

function test_nominalScenarioConstructor
    scenario = matRad_NominalScenario();
    assertTrue(isa(scenario, 'matRad_NominalScenario'));
    assertTrue(isa(scenario, 'matRad_ScenarioModel'));
    assertEqual(scenario.name, 'nomScen');
    assertEqual(scenario.ctScenProb, [1 1]);
    assertEqual(scenario.numOfCtScen, 1);
    assertEqual(scenario.totNumScen, 1);
    assertEqual(scenario.totNumShiftScen, 1);
    assertEqual(scenario.totNumRangeScen, 1);
    assertEqual(scenario.relRangeShift, 0);
    assertEqual(scenario.absRangeShift, 0);
    assertEqual(scenario.isoShift, zeros(1,3));
    assertEqual(scenario.maxAbsRangeShift, 0);
    assertEqual(scenario.maxRelRangeShift, 0);
    assertEqual(scenario.scenMask, true(1,1,1));
    assertEqual(scenario.linearMask, [1 1 1]);
    assertElementsAlmostEqual(scenario.scenProb,helper_mvarGauss(scenario));
    assertEqual(scenario.scenForProb,[1 zeros(1,5)]);
    assertEqual(scenario.scenWeight, 1);

function test_nominalScenarioConstructorWithCt
    n = 5;
    ct = struct('numOfCtScen',n);
    scenario = matRad_NominalScenario(ct);
    assertTrue(isa(scenario, 'matRad_NominalScenario'));
    assertTrue(isa(scenario, 'matRad_ScenarioModel'));
    assertEqual(scenario.name, 'nomScen');
    assertEqual(scenario.ctScenProb, [(1:n)' ones(n,1)./n]);
    assertEqual(scenario.numOfCtScen, n);
    assertEqual(scenario.totNumScen, n);
    assertEqual(scenario.totNumShiftScen, 1);
    assertEqual(scenario.totNumRangeScen, 1);
    assertEqual(scenario.relRangeShift, zeros(n,1));
    assertEqual(scenario.absRangeShift, zeros(n,1));
    assertEqual(scenario.isoShift, zeros(n,3));
    assertEqual(scenario.maxAbsRangeShift, 0);
    assertEqual(scenario.maxRelRangeShift, 0);
    assertTrue(isequal(scenario.scenMask, true(n,1,1)));
    assertEqual(scenario.linearMask, [(1:n)' ones(n,2)]);
    assertElementsAlmostEqual(scenario.scenProb,helper_mvarGauss(scenario));
    assertEqual(scenario.scenForProb,[(1:n)' zeros(n,5)]);
    assertEqual(scenario.scenWeight, ones(n,1)./n);

function test_nominalScenarioExtractSingleScenario
    scenario = matRad_NominalScenario();
    newInstance = scenario.extractSingleScenario(1);
    assertTrue(isequal(scenario,newInstance));

function test_nominalScenarioExtractSingleScenarioWithCtScen
    n = 5;
    ct = struct('numOfCtScen',n);
    refScen = matRad_NominalScenario(ct);

    for scenNum = 1:refScen.totNumScen
        scenario = refScen.extractSingleScenario(scenNum);
        assertTrue(isa(scenario, 'matRad_NominalScenario'));
        ctScenIx = refScen.ctScenIx(scenNum);
        ctScenNum = find(ctScenIx == refScen.ctScenProb(:,1));
        assertEqual(scenario.ctScenProb, refScen.ctScenProb(ctScenNum,:));
        assertEqual(scenario.numOfCtScen, 1);
        assertEqual(scenario.totNumScen, 1);
        assertEqual(scenario.totNumShiftScen, 1);
        assertEqual(scenario.totNumRangeScen, 1);
        assertEqual(scenario.relRangeShift, refScen.relRangeShift(scenNum));
        assertEqual(scenario.absRangeShift, refScen.absRangeShift(scenNum));
        assertEqual(scenario.isoShift, refScen.isoShift(scenNum,:));
        assertEqual(scenario.maxAbsRangeShift, max(abs(refScen.absRangeShift(scenNum))));
        assertEqual(scenario.maxRelRangeShift, max(abs(refScen.relRangeShift(scenNum))));
        assertTrue(scenario.scenMask(ctScenIx,1,1));
        assertTrue(numel(find(scenario.scenMask)) == 1);
        assertEqual(scenario.linearMask, [ctScenIx 1 1]);
        assertElementsAlmostEqual(scenario.scenProb,helper_mvarGauss(scenario));
        assertEqual(scenario.scenForProb,refScen.scenForProb(scenNum,:));
        assertEqual(scenario.scenWeight, refScen.scenWeight(scenNum));
    end


