function test_suite = test_randomScenarios

test_functions=localfunctions();

initTestSuite;

function test_randomScenarioConstructor
    scenario = matRad_RandomScenarios();

    %Defaults
    nSamples = 10;

    assertTrue(isa(scenario, 'matRad_RandomScenarios'));
    assertTrue(isa(scenario, 'matRad_ScenarioModel'));
    assertEqual(scenario.name, 'rndScen');
    %Test correct standard values & sizes
    assertEqual(scenario.ctScenProb, [1 1]);
    assertEqual(scenario.numOfCtScen, 1);
    assertEqual(scenario.totNumScen, nSamples); 
    assertEqual(scenario.totNumShiftScen, nSamples);
    assertEqual(scenario.totNumRangeScen, nSamples);
    assertEqual(size(scenario.relRangeShift),[scenario.totNumScen,1]);
    assertEqual(size(scenario.absRangeShift),[scenario.totNumScen,1]);
    assertEqual(size(scenario.isoShift),[scenario.totNumScen,3]);
    %assertEqual(scenario.isoShift, 2.25 * ones(1,3));
    assertEqual(scenario.maxAbsRangeShift, max(scenario.absRangeShift));
    assertEqual(scenario.maxRelRangeShift, max(scenario.relRangeShift));
    assertEqual(size(scenario.scenMask), [scenario.numOfCtScen,scenario.totNumShiftScen,scenario.totNumRangeScen]);
    %assertEqual(scenario.scenMask, true(1,1,1));
    assertEqual(size(scenario.linearMask), [scenario.totNumScen,3]);

    tmpScenMask = permute(scenario.scenMask,[2 3 1]);
    [tmp(:,2),tmp(:,3),tmp(:,1)] = ind2sub(size(tmpScenMask),find(tmpScenMask));
    assertEqual(tmp,scenario.linearMask);

    %assertEqual(ind2sub(find()))
    %assertEqual(scenario.linearMask, [1 1 1]);
    assertElementsAlmostEqual(scenario.scenProb,helper_mvarGauss(scenario));

    tmp = [scenario.ctScenIx scenario.isoShift scenario.absRangeShift scenario.relRangeShift];    
    assertEqual(scenario.scenForProb,tmp);

    assertEqual(numel(unique(scenario.scenWeight)),scenario.numOfCtScen);    

function test_randomScenarioConstructorWithCt
    n = 5;
    ct = struct('numOfCtScen',n);

    %Defaults
    nSamples = 10;

    scenario = matRad_RandomScenarios(ct);
    assertTrue(isa(scenario, 'matRad_RandomScenarios'));
    assertTrue(isa(scenario, 'matRad_ScenarioModel'));
    assertEqual(scenario.name, 'rndScen');
    %Test correct standard values & sizes
    assertEqual(scenario.ctScenProb, [(1:n)' ones(n,1)./n]);
    assertEqual(scenario.numOfCtScen, n);
    assertEqual(scenario.totNumScen, nSamples*n);
    assertEqual(scenario.totNumShiftScen, nSamples);
    assertEqual(scenario.totNumRangeScen, nSamples);
    assertEqual(size(scenario.relRangeShift),[scenario.totNumScen,1]);
    assertEqual(size(scenario.absRangeShift),[scenario.totNumScen,1]);
    assertEqual(size(scenario.isoShift),[scenario.totNumScen,3]);
    %assertEqual(scenario.isoShift, 2.25 * ones(1,3));
    assertEqual(scenario.maxAbsRangeShift, max(scenario.absRangeShift));
    assertEqual(scenario.maxRelRangeShift, max(scenario.relRangeShift));
    assertEqual(size(scenario.scenMask), [scenario.numOfCtScen,scenario.totNumShiftScen,scenario.totNumRangeScen]);
    %assertEqual(scenario.scenMask, true(1,1,1));
    assertEqual(size(scenario.linearMask), [scenario.totNumScen,3]);
    
    tmpScenMask = permute(scenario.scenMask,[2 3 1]);
    [tmp(:,2),tmp(:,3),tmp(:,1)] = ind2sub(size(tmpScenMask),find(tmpScenMask));
    assertEqual(tmp,scenario.linearMask);

    %assertEqual(ind2sub(find()))
    %assertEqual(scenario.linearMask, [1 1 1]);
    assertElementsAlmostEqual(scenario.scenProb,helper_mvarGauss(scenario));

    tmp = [scenario.ctScenIx scenario.isoShift scenario.absRangeShift scenario.relRangeShift];    
    assertEqual(scenario.scenForProb,tmp);
    assertEqual(numel(unique(scenario.scenWeight)),numel(unique(scenario.ctScenProb(:,2))));   

    
function test_randomScenarioExtractSingleScenario
    refScen = matRad_RandomScenarios();
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


function test_randomScenarioExtractSingleScenarioWithCtScen
    n = 5;
    ct = struct('numOfCtScen',n);
    refScen = matRad_RandomScenarios(ct);
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

