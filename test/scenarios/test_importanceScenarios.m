function test_suite = test_importanceScenarios

test_functions=localfunctions();

initTestSuite;

function test_importanceScenarioConstructor
    scenario = matRad_ImportanceScenarios();

    %Defaults
    nShift = 25;
    nRange = 9;
    nTot = nShift + nRange - 1; %Does include the nominal scenario by default 


    assertTrue(isa(scenario, 'matRad_ImportanceScenarios'));
    assertTrue(isa(scenario, 'matRad_GriddedScenariosAbstract'));
    assertTrue(isa(scenario, 'matRad_ScenarioModel'));
    assertEqual(scenario.name, 'impScen');
    %Test correct standard values & sizes
    assertEqual(scenario.ctScenProb, [1 1]);
    assertEqual(scenario.numOfCtScen, 1);
    assertEqual(scenario.totNumScen, nTot); 
    assertEqual(scenario.totNumShiftScen, nShift);
    assertEqual(scenario.totNumRangeScen, nRange);
    assertEqual(size(scenario.relRangeShift),[scenario.totNumScen,1]);
    assertEqual(size(scenario.absRangeShift),[scenario.totNumScen,1]);
    assertEqual(size(scenario.isoShift),[scenario.totNumScen,3]);
    %assertEqual(scenario.isoShift, 2.25 * ones(1,3));
    assertEqual(scenario.maxAbsRangeShift, max(scenario.absRangeShift));
    assertEqual(scenario.maxRelRangeShift, max(scenario.relRangeShift));
    assertEqual(size(scenario.scenMask), [scenario.numOfCtScen,scenario.totNumShiftScen,scenario.totNumRangeScen]);
    %assertEqual(scenario.scenMask, true(1,1,1));
    assertEqual(size(scenario.linearMask), [scenario.totNumScen,3]);

    [tmp(:,1),tmp(:,2),tmp(:,3)] = ind2sub(size(scenario.scenMask),find(scenario.scenMask));
    assertEqual(tmp,scenario.linearMask);

    %assertEqual(ind2sub(find()))
    %assertEqual(scenario.linearMask, [1 1 1]);
    assertElementsAlmostEqual(scenario.scenProb,helper_mvarGauss(scenario));

    tmp = [scenario.ctScenIx scenario.isoShift scenario.absRangeShift scenario.relRangeShift];    
    assertEqual(scenario.scenForProb,tmp);
    assertEqual(scenario.scenWeight,scenario.scenProb./sum(scenario.scenProb));

function test_importanceScenarioConstructorWithCt
    n = 5;
    ct = struct('numOfCtScen',n);

    %Defaults
    nShift = 25;
    nRange = 9;
    nTot = nShift + nRange - 1; %Does include the nominal scenario by default 

    scenario = matRad_ImportanceScenarios(ct);
    assertTrue(isa(scenario, 'matRad_ImportanceScenarios'));
    assertTrue(isa(scenario, 'matRad_GriddedScenariosAbstract'));
    assertTrue(isa(scenario, 'matRad_ScenarioModel'));
    assertEqual(scenario.name, 'impScen');
    %Test correct standard values & sizes
    assertEqual(scenario.ctScenProb, [(1:n)' ones(n,1)./n]);
    assertEqual(scenario.numOfCtScen, n);
    assertEqual(scenario.totNumScen, nTot*n);
    assertEqual(scenario.totNumShiftScen, nShift);
    assertEqual(scenario.totNumRangeScen, nRange);
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
    assertEqual(scenario.scenWeight,scenario.scenProb./sum(scenario.scenProb));

    
function test_importanceScenarioExtractSingleScenario
    refScen = matRad_ImportanceScenarios();
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


function test_importanceScenarioExtractSingleScenarioWithCtScen
    n = 5;
    ct = struct('numOfCtScen',n);
    refScen = matRad_ImportanceScenarios(ct);

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
    

    function test_importanceScenarioCombineRange

        model = matRad_ImportanceScenarios();
    
        assertExceptionThrown(@() helper_assignmentTest(model,'combineRange','hello'),'matRad:Error');
        assertTrue(model.combineRange);
        
        nRangeScen = model.totNumRangeScen;

        model.combineRange = false;
        assertFalse(model.combineRange);
        assertEqual(model.totNumRangeScen,nRangeScen^2);
        assertEqual(model.totNumScen,model.totNumRangeScen - 1 + model.totNumShiftScen);
    
    function test_importanceScenarioShiftCombinations
    
        model = matRad_ImportanceScenarios();
    
        assertExceptionThrown(@() helper_assignmentTest(model,'combinations','hello'),'matRad:Error');
        assertEqual(model.combinations,'none');
    
        nSetupPoints = model.numOfSetupGridPoints;
        nRangePoints = model.numOfRangeGridPoints;

        model.combinations = 'shift';
        assertEqual(model.combinations,'shift');
        assertEqual(model.totNumShiftScen,nSetupPoints^3);
        assertEqual(model.totNumScen,model.totNumRangeScen - 1 + model.totNumShiftScen);
    
        model.combinations = 'all';
        assertEqual(model.combinations,'all');
        assertEqual(model.totNumShiftScen,nSetupPoints^3);
        assertEqual(model.totNumRangeScen,nRangePoints);
        assertEqual(model.totNumScen,model.totNumRangeScen * model.totNumShiftScen);
    
        model.combineRange = false;
        assertEqual(model.totNumShiftScen,nSetupPoints^3);
        assertEqual(model.totNumRangeScen,nRangePoints^2);
        assertEqual(model.totNumScen,model.totNumRangeScen * model.totNumShiftScen);
    
        model.combinations = 'shift';
        assertEqual(model.totNumShiftScen,nSetupPoints^3);
        assertEqual(model.totNumRangeScen,nRangePoints^2);
        assertEqual(model.totNumScen,model.totNumRangeScen + model.totNumShiftScen - 1);