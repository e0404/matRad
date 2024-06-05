function test_suite = test_nominalScenario

test_functions=localfunctions();

initTestSuite;

function p = helper_mvarGauss(model)
    Sigma = diag([model.shiftSD,model.rangeAbsSD,model.rangeRelSD./100].^2);
    d = size(Sigma,1);
    [cs,~] = chol(Sigma);
    
    % Compute from Gaussian errors
    p = (2*pi)^(-d/2) * exp(-0.5*sum((model.scenForProb(:,2:end)/cs).^2, 2)) / prod(diag(cs));

    % Now multiplay with the phase probability
    tmpPhaseProb = arrayfun(@(phase) model.phaseProbability(phase),model.scenForProb(:,1));
    p = p .* tmpPhaseProb;


function test_nominalScenarioConstructor
    scenario = matRad_NominalScenario();
    assertTrue(isa(scenario, 'matRad_NominalScenario'));
    assertTrue(isa(scenario, 'matRad_ScenarioModel'));
    assertEqual(scenario.name, 'nomScen');
    assertEqual(scenario.phaseProbability, 1);
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
    assertEqual(scenario.phaseProbability, ones(n,1)./n);
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
    scenNum = 1;
    ct = struct('numOfCtScen',n);
    refScen = matRad_NominalScenario(ct);
    scenario = refScen.extractSingleScenario(scenNum);

    assertTrue(isa(scenario, 'matRad_NominalScenario'));
    assertTrue(isa(scenario, 'matRad_ScenarioModel'));
    assertEqual(scenario.name, 'nomScen');
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


