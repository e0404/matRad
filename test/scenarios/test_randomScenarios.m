function test_suite = test_randomScenarios

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

function test_randomScenarioConstructor
    scenario = matRad_RandomScenarios();

    %Defaults
    nSamples = 10;

    assertTrue(isa(scenario, 'matRad_RandomScenarios'));
    assertTrue(isa(scenario, 'matRad_ScenarioModel'));
    assertEqual(scenario.name, 'rndScen');
    %Test correct standard values & sizes
    assertEqual(scenario.phaseProbability, 1);
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

    tmp = [scenario.ctScen scenario.isoShift scenario.absRangeShift scenario.relRangeShift];    
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
    assertEqual(scenario.phaseProbability, ones(n,1)./n);
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

    tmp = [scenario.ctScen scenario.isoShift scenario.absRangeShift scenario.relRangeShift];    
    assertEqual(scenario.scenForProb,tmp);
    assertEqual(numel(unique(scenario.scenWeight)),numel(unique(scenario.phaseProbability)));   

    
function test_randomScenarioExtractSingleScenario
    scenNum = 1;
    refScen = matRad_RandomScenarios();
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


function test_randomScenarioExtractSingleScenarioWithCtScen
    n = 5;
    scenNum = 1;
    ct = struct('numOfCtScen',n);
    refScen = matRad_RandomScenarios(ct);
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
    

