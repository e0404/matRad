function test_suite = test_wcScenarios

test_functions=localfunctions();

initTestSuite;

function assignmentTestHelper(model,property,value)
    model.(property) = value;

function p = helper_mvarGauss(model)
    Sigma = diag([model.shiftSD,model.rangeAbsSD,model.rangeRelSD./100].^2);
    d = size(Sigma,1);
    [cs,~] = chol(Sigma);
    
    % Compute from Gaussian errors
    p = (2*pi)^(-d/2) * exp(-0.5*sum((model.scenForProb(:,2:end)/cs).^2, 2)) / prod(diag(cs));

    % Now multiplay with the phase probability
    tmpPhaseProb = arrayfun(@(phase) model.phaseProbability(phase),model.scenForProb(:,1));
    p = p .* tmpPhaseProb;
    


function test_worstCaseScenarioConstructor
    scenario = matRad_WorstCaseScenarios();
    assertTrue(isa(scenario, 'matRad_WorstCaseScenarios'));
    assertTrue(isa(scenario, 'matRad_GriddedScenariosAbstract'));
    assertTrue(isa(scenario, 'matRad_ScenarioModel'));
    assertEqual(scenario.name, 'wcScen');
    %Test correct standard values & sizes
    assertEqual(scenario.phaseProbability, 1);
    assertEqual(scenario.numOfCtScen, 1);
    assertEqual(scenario.totNumScen, 9);
    assertEqual(scenario.totNumShiftScen, 7);
    assertEqual(scenario.totNumRangeScen, 3);
    assertEqual(size(scenario.relRangeShift),[scenario.totNumScen,1]);
    assertEqual(size(scenario.absRangeShift),[scenario.totNumScen,1]);
    assertEqual(size(scenario.isoShift),[scenario.totNumScen,3]);
    assertEqual(scenario.relRangeShift, [zeros(7,1); -0.035; 0.035]);
    assertEqual(scenario.absRangeShift, [zeros(7,1); -1; 1]);
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

    tmp = [scenario.ctScen scenario.isoShift scenario.absRangeShift scenario.relRangeShift];    
    assertEqual(scenario.scenForProb,tmp);
    assertEqual(scenario.scenWeight,scenario.scenProb./sum(scenario.scenProb));

function test_worstCaseScenarioConstructorWithCt
    n = 5;
    ct = struct('numOfCtScen',n);

    scenario = matRad_WorstCaseScenarios(ct);
    assertTrue(isa(scenario, 'matRad_WorstCaseScenarios'));
    assertTrue(isa(scenario, 'matRad_GriddedScenariosAbstract'));
    assertTrue(isa(scenario, 'matRad_ScenarioModel'));
    assertEqual(scenario.name, 'wcScen');
    %Test correct standard values & sizes
    assertEqual(scenario.phaseProbability, ones(n,1)./n);
    assertEqual(scenario.numOfCtScen, n);
    assertEqual(scenario.totNumScen, 9*n);
    assertEqual(scenario.totNumShiftScen, 7);
    assertEqual(scenario.totNumRangeScen, 3);
    assertEqual(size(scenario.relRangeShift),[scenario.totNumScen,1]);
    assertEqual(size(scenario.absRangeShift),[scenario.totNumScen,1]);
    assertEqual(size(scenario.isoShift),[scenario.totNumScen,3]);
    assertEqual(scenario.relRangeShift, repmat([zeros(7,1); -0.035; 0.035],n,1));
    assertEqual(scenario.absRangeShift, repmat([zeros(7,1); -1; 1],n,1));
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

    tmp = [scenario.ctScen scenario.isoShift scenario.absRangeShift scenario.relRangeShift];    
    assertEqual(scenario.scenForProb,tmp);
    assertEqual(scenario.scenWeight,scenario.scenProb./sum(scenario.scenProb));

    
function test_worstCaseScenarioExtractSingleScenario
    scenNum = 1;
    refScen = matRad_WorstCaseScenarios();
    scenario = refScen.extractSingleScenario(scenNum);
    assertTrue(isa(scenario, 'matRad_NominalScenario'));
    assertEqual(scenario.phaseProbability, refScen.phaseProbability(scenNum));
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
    assertEqual(scenario.scenWeight, refScen.scenWeight(scenNum));


function test_worstCaseScenarioExtractSingleScenarioWithCtScen
    n = 5;
    scenNum = 1;
    ct = struct('numOfCtScen',n);
    refScen = matRad_WorstCaseScenarios(ct);
    scenario = refScen.extractSingleScenario(1);

    assertTrue(isa(scenario, 'matRad_NominalScenario'));
    assertEqual(scenario.phaseProbability, refScen.phaseProbability(scenNum));
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
    assertEqual(scenario.scenWeight, refScen.scenWeight(scenNum));
    
function test_worstCaseScenarioCombineRange

    model = matRad_WorstCaseScenarios();

    assertExceptionThrown(@() assignmentTestHelper(model,'combineRange','hello'),'matRad:Error');
    assertTrue(model.combineRange);

    assertEqual(model.totNumRangeScen,3);
    model.combineRange = false;
    assertFalse(model.combineRange);
    assertEqual(model.totNumScen,15);
    assertEqual(model.totNumRangeScen,9);

function test_worstCaseScenarioShiftCombinations

    model = matRad_WorstCaseScenarios();

    assertExceptionThrown(@() assignmentTestHelper(model,'combinations','hello'),'matRad:Error');
    assertEqual(model.combinations,'none');

    model.combinations = 'shift';
    assertEqual(model.combinations,'shift');
    assertEqual(model.totNumShiftScen,27);
    assertEqual(model.totNumRangeScen,3);
    assertEqual(model.totNumScen,29);

    model.combinations = 'all';
    assertEqual(model.combinations,'all');
    assertEqual(model.totNumShiftScen,27);
    assertEqual(model.totNumRangeScen,3);
    assertEqual(model.totNumScen,81);

    model.combineRange = false;
    assertEqual(model.totNumShiftScen,27);
    assertEqual(model.totNumRangeScen,9);
    assertEqual(model.totNumScen,243);

    model.combinations = 'shift';
    assertEqual(model.totNumShiftScen,27);
    assertEqual(model.totNumRangeScen,9);
    assertEqual(model.totNumScen,35);