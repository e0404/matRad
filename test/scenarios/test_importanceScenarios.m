function test_suite = test_importanceScenarios

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
    assertEqual(scenario.phaseProbability, 1);
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

    tmp = [scenario.ctScen scenario.isoShift scenario.absRangeShift scenario.relRangeShift];    
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
    assertEqual(scenario.phaseProbability, ones(n,1)./n);
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

    tmp = [scenario.ctScen scenario.isoShift scenario.absRangeShift scenario.relRangeShift];    
    assertEqual(scenario.scenForProb,tmp);
    assertEqual(scenario.scenWeight,scenario.scenProb./sum(scenario.scenProb));

    
function test_importanceScenarioExtractSingleScenario
    refScen = matRad_ImportanceScenarios();
    scenario = refScen.extractSingleScenario(1);
    assertTrue(isa(scenario, 'matRad_NominalScenario'));
    assertEqual(scenario.phaseProbability, refScen.phaseProbability(1));
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


function test_importanceScenarioExtractSingleScenarioWithCtScen
    n = 5;
    ct = struct('numOfCtScen',n);
    refScen = matRad_ImportanceScenarios(ct);
    scenario = refScen.extractSingleScenario(1);

    assertTrue(isa(scenario, 'matRad_NominalScenario'));
    assertEqual(scenario.phaseProbability, refScen.phaseProbability(1));
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
    

    function test_importanceScenarioCombineRange

        model = matRad_ImportanceScenarios();
    
        assertExceptionThrown(@() assignmentTestHelper(model,'combineRange','hello'),'matRad:Error');
        assertTrue(model.combineRange);
        
        nRangeScen = model.totNumRangeScen;

        model.combineRange = false;
        assertFalse(model.combineRange);
        assertEqual(model.totNumRangeScen,nRangeScen^2);
        assertEqual(model.totNumScen,model.totNumRangeScen - 1 + model.totNumShiftScen);
    
    function test_importanceScenarioShiftCombinations
    
        model = matRad_ImportanceScenarios();
    
        assertExceptionThrown(@() assignmentTestHelper(model,'combinations','hello'),'matRad:Error');
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