function test_suite = test_wcScenarios

test_functions = localfunctions();

initTestSuite;

function test_worstCaseScenarioConstructor
scenario = matRad_WorstCaseScenarios();
assertTrue(isa(scenario, 'matRad_WorstCaseScenarios'));
assertTrue(isa(scenario, 'matRad_GriddedScenariosAbstract'));
assertTrue(isa(scenario, 'matRad_ScenarioModel'));
assertEqual(scenario.shortName, 'wcScen');
% Test correct standard values & sizes
assertEqual(scenario.ctScenProb, [1 1]);
assertEqual(scenario.numOfCtScen, 1);
assertEqual(scenario.totNumScen, 9);
assertEqual(scenario.totNumShiftScen, 7);
assertEqual(scenario.totNumRangeScen, 3);
assertEqual(size(scenario.relRangeShift), [scenario.totNumScen, 1]);
assertEqual(size(scenario.absRangeShift), [scenario.totNumScen, 1]);
assertEqual(size(scenario.isoShift), [scenario.totNumScen, 3]);
assertEqual(scenario.relRangeShift, [zeros(7, 1); -0.035; 0.035]);
assertEqual(scenario.absRangeShift, [zeros(7, 1); -1; 1]);
% assertEqual(scenario.isoShift, 2.25 * ones(1,3));
assertEqual(scenario.maxAbsRangeShift, max(scenario.absRangeShift));
assertEqual(scenario.maxRelRangeShift, max(scenario.relRangeShift));
assertEqual(size(scenario.scenMask), [scenario.numOfCtScen, scenario.totNumShiftScen, scenario.totNumRangeScen]);
% assertEqual(scenario.scenMask, true(1,1,1));
assertEqual(size(scenario.linearMask), [scenario.totNumScen, 3]);

[tmp(:, 1), tmp(:, 2), tmp(:, 3)] = ind2sub(size(scenario.scenMask), find(scenario.scenMask));
assertEqual(tmp, scenario.linearMask);

% assertEqual(ind2sub(find()))
% assertEqual(scenario.linearMask, [1 1 1]);
assertElementsAlmostEqual(scenario.scenProb, helper_mvarGauss(scenario));

tmp = [scenario.ctScenIx scenario.isoShift scenario.absRangeShift scenario.relRangeShift];
assertEqual(scenario.scenForProb, tmp);
assertEqual(scenario.scenWeight, scenario.scenProb ./ sum(scenario.scenProb));

function test_worstCaseScenarioConstructorWithCt
n = 5;
ct = struct('numOfCtScen', n);

scenario = matRad_WorstCaseScenarios(ct);
assertTrue(isa(scenario, 'matRad_WorstCaseScenarios'));
assertTrue(isa(scenario, 'matRad_GriddedScenariosAbstract'));
assertTrue(isa(scenario, 'matRad_ScenarioModel'));
assertEqual(scenario.shortName, 'wcScen');
% Test correct standard values & sizes
assertEqual(scenario.ctScenProb, [(1:n)' ones(n, 1) ./ n]);
assertEqual(scenario.numOfCtScen, n);
assertEqual(scenario.totNumScen, 9 * n);
assertEqual(scenario.totNumShiftScen, 7);
assertEqual(scenario.totNumRangeScen, 3);
assertEqual(size(scenario.relRangeShift), [scenario.totNumScen, 1]);
assertEqual(size(scenario.absRangeShift), [scenario.totNumScen, 1]);
assertEqual(size(scenario.isoShift), [scenario.totNumScen, 3]);
assertEqual(scenario.relRangeShift, repmat([zeros(7, 1); -0.035; 0.035], n, 1));
assertEqual(scenario.absRangeShift, repmat([zeros(7, 1); -1; 1], n, 1));
% assertEqual(scenario.isoShift, 2.25 * ones(1,3));
assertEqual(scenario.maxAbsRangeShift, max(scenario.absRangeShift));
assertEqual(scenario.maxRelRangeShift, max(scenario.relRangeShift));
assertEqual(size(scenario.scenMask), [scenario.numOfCtScen, scenario.totNumShiftScen, scenario.totNumRangeScen]);
% assertEqual(scenario.scenMask, true(1,1,1));
assertEqual(size(scenario.linearMask), [scenario.totNumScen, 3]);

tmpScenMask = permute(scenario.scenMask, [2 3 1]);

[tmp(:, 2), tmp(:, 3), tmp(:, 1)] = ind2sub(size(tmpScenMask), find(tmpScenMask));
assertEqual(tmp, scenario.linearMask);

% assertEqual(ind2sub(find()))
% assertEqual(scenario.linearMask, [1 1 1]);
assertElementsAlmostEqual(scenario.scenProb, helper_mvarGauss(scenario));

tmp = [scenario.ctScenIx scenario.isoShift scenario.absRangeShift scenario.relRangeShift];
assertEqual(scenario.scenForProb, tmp);
assertEqual(scenario.scenWeight, scenario.scenProb ./ sum(scenario.scenProb));

function test_worstCaseScenarioExtractSingleScenario
refScen = matRad_WorstCaseScenarios();

for scenNum = 1:refScen.totNumScen
    scenario = refScen.extractSingleScenario(scenNum);
    assertTrue(isa(scenario, 'matRad_NominalScenario'));
    ctScenIx = refScen.ctScenIx(scenNum);
    ctScenNum = find(ctScenIx == refScen.ctScenProb(:, 1));
    assertEqual(scenario.ctScenProb, refScen.ctScenProb(ctScenNum, :));
    assertEqual(scenario.numOfCtScen, 1);
    assertEqual(scenario.totNumScen, 1);
    assertEqual(scenario.totNumShiftScen, 1);
    assertEqual(scenario.totNumRangeScen, 1);
    assertEqual(scenario.relRangeShift, refScen.relRangeShift(scenNum));
    assertEqual(scenario.absRangeShift, refScen.absRangeShift(scenNum));
    assertEqual(scenario.isoShift, refScen.isoShift(scenNum, :));
    assertEqual(scenario.maxAbsRangeShift, max(abs(refScen.absRangeShift(scenNum))));
    assertEqual(scenario.maxRelRangeShift, max(abs(refScen.relRangeShift(scenNum))));
    assertTrue(scenario.scenMask(ctScenIx, 1, 1));
    assertTrue(numel(find(scenario.scenMask)) == 1);
    assertEqual(scenario.linearMask, [ctScenIx 1 1]);
    assertElementsAlmostEqual(scenario.scenProb, helper_mvarGauss(scenario));
    assertEqual(scenario.scenForProb, refScen.scenForProb(scenNum, :));
    assertEqual(scenario.scenWeight, refScen.scenWeight(scenNum));
end

function test_worstCaseScenarioExtractSingleScenarioWithCtScen
n = 5;
scenNum = 1;
ct = struct('numOfCtScen', n);
refScen = matRad_WorstCaseScenarios(ct);
for scenNum = 1:refScen.totNumScen
    scenario = refScen.extractSingleScenario(scenNum);
    assertTrue(isa(scenario, 'matRad_NominalScenario'));
    ctScenIx = refScen.ctScenIx(scenNum);
    ctScenNum = find(ctScenIx == refScen.ctScenProb(:, 1));
    assertEqual(scenario.ctScenProb, refScen.ctScenProb(ctScenNum, :));
    assertEqual(scenario.numOfCtScen, 1);
    assertEqual(scenario.totNumScen, 1);
    assertEqual(scenario.totNumShiftScen, 1);
    assertEqual(scenario.totNumRangeScen, 1);
    assertEqual(scenario.relRangeShift, refScen.relRangeShift(scenNum));
    assertEqual(scenario.absRangeShift, refScen.absRangeShift(scenNum));
    assertEqual(scenario.isoShift, refScen.isoShift(scenNum, :));
    assertEqual(scenario.maxAbsRangeShift, max(abs(refScen.absRangeShift(scenNum))));
    assertEqual(scenario.maxRelRangeShift, max(abs(refScen.relRangeShift(scenNum))));
    assertTrue(scenario.scenMask(ctScenIx, 1, 1));
    assertTrue(numel(find(scenario.scenMask)) == 1);
    assertEqual(scenario.linearMask, [ctScenIx 1 1]);
    assertElementsAlmostEqual(scenario.scenProb, helper_mvarGauss(scenario));
    assertEqual(scenario.scenForProb, refScen.scenForProb(scenNum, :));
    assertEqual(scenario.scenWeight, refScen.scenWeight(scenNum));
end

function test_worstCaseScenarioCombineRange

model = matRad_WorstCaseScenarios();

assertExceptionThrown(@() helper_assignmentTest(model, 'combineRange', 'hello'), 'matRad:Error');
assertTrue(model.combineRange);

assertEqual(model.totNumRangeScen, 3);
model.combineRange = false;
assertFalse(model.combineRange);
assertEqual(model.totNumScen, 15);
assertEqual(model.totNumRangeScen, 9);

function test_worstCaseScenarioShiftCombinations

model = matRad_WorstCaseScenarios();

assertExceptionThrown(@() helper_assignmentTest(model, 'combinations', 'hello'), 'matRad:Error');
assertEqual(model.combinations, 'none');

model.combinations = 'shift';
assertEqual(model.combinations, 'shift');
assertEqual(model.totNumShiftScen, 27);
assertEqual(model.totNumRangeScen, 3);
assertEqual(model.totNumScen, 29);

model.combinations = 'all';
assertEqual(model.combinations, 'all');
assertEqual(model.totNumShiftScen, 27);
assertEqual(model.totNumRangeScen, 3);
assertEqual(model.totNumScen, 81);

model.combineRange = false;
assertEqual(model.totNumShiftScen, 27);
assertEqual(model.totNumRangeScen, 9);
assertEqual(model.totNumScen, 243);

model.combinations = 'shift';
assertEqual(model.totNumShiftScen, 27);
assertEqual(model.totNumRangeScen, 9);
assertEqual(model.totNumScen, 35);

function test_worstCaseScenarioGridPoints
% Default: numOfSetupGridPoints=3, numOfRangeGridPoints=3 -> 9 scenarios
% (1 nominal + 6 shifts = 7 shift scenarios, 3 range scenarios)
model = matRad_WorstCaseScenarios();
assertEqual(model.numOfSetupGridPoints, 3);
assertEqual(model.numOfRangeGridPoints, 3);
assertEqual(model.totNumShiftScen, 7);
assertEqual(model.totNumRangeScen, 3);
assertEqual(model.totNumScen, 9);

% Invalid grid point values should throw an error
assertExceptionThrown(@() helper_assignmentTest(model, 'numOfSetupGridPoints', 2), 'matRad:Error');
assertExceptionThrown(@() helper_assignmentTest(model, 'numOfRangeGridPoints', 2), 'matRad:Error');

% numOfRangeGridPoints=1: range worst cases removed -> 7 scenarios
% (7 shift scenarios, 1 range scenario)
model.numOfRangeGridPoints = 1;
assertEqual(model.numOfRangeGridPoints, 1);
assertEqual(model.totNumShiftScen, 7);
assertEqual(model.totNumRangeScen, 1);
assertEqual(model.totNumScen, 7);

% numOfSetupGridPoints=1: shift worst cases removed -> 3 scenarios
% (1 shift scenario, 1 range scenario -> but range is still 1 from above)
% Reset range first, then set setup to 1
model.numOfRangeGridPoints = 3;
model.numOfSetupGridPoints = 1;
assertEqual(model.numOfSetupGridPoints, 1);
assertEqual(model.totNumShiftScen, 1);
assertEqual(model.totNumRangeScen, 3);
assertEqual(model.totNumScen, 3);

function test_worstCaseScenarioSub2scenIxDegenerate
% Regression test for sub2scenIx when scenMask dimensions are degenerate
% (i.e., MATLAB squeezes trailing/leading singletons).

% Case 1: numOfRangeGridPoints=1 -> scenMask is [1 x N_S x 1], MATLAB
% squeezes to [1 x N_S] (row vector). sub2scenIx must NOT fall into the
% CT-only branch (iscolumn = false for row vectors).
model = matRad_WorstCaseScenarios();
model.numOfRangeGridPoints = 1;
assertEqual(model.totNumRangeScen, 1);
assertEqual(model.totNumShiftScen, 7);
assertEqual(model.numOfCtScen, 1);
% Verify sub2scenIx round-trips correctly for every scenario
for s = 1:model.totNumScen
    scenIx = model.sub2scenIx(model.linearMask(s, 1), model.linearMask(s, 2), model.linearMask(s, 3));
    assertEqual(model.scenNum(scenIx), s);
end

% Case 2: numOfSetupGridPoints=1 -> scenMask is [1 x 1 x N_R], stays 3-D.
% sub2ind path must work with the known dimensions, not size(scenMask).
model2 = matRad_WorstCaseScenarios();
model2.numOfSetupGridPoints = 1;
assertEqual(model2.totNumShiftScen, 1);
assertEqual(model2.totNumRangeScen, 3);
for s = 1:model2.totNumScen
    scenIx = model2.sub2scenIx(model2.linearMask(s, 1), model2.linearMask(s, 2), model2.linearMask(s, 3));
    assertEqual(model2.scenNum(scenIx), s);
end
