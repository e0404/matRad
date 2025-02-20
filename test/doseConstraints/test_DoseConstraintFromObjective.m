function test_suite = test_DoseConstraintFromObjective

test_functions=localfunctions();

initTestSuite;

function test_matRad_DoseConstraintFromObjective
    objective = DoseObjectives.matRad_SquaredOverdosing(100,30);
    obj = DoseConstraints.matRad_DoseConstraintFromObjective(objective);   %default
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseConstraints.matRad_DoseConstraintFromObjective'));
    assertEqual(obj.parameters{1},1e-5);
    assertEqual(obj.parameters{2},1e-3);
    
    
    obj = DoseConstraints.matRad_DoseConstraintFromObjective(objective,1e-3,1e-3);
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseConstraints.matRad_DoseConstraintFromObjective'));

function test_matRad_DoseConstraintFromObjective_bounds
    objective = DoseObjectives.matRad_SquaredOverdosing(100,30);
    obj = DoseConstraints.matRad_DoseConstraintFromObjective(objective,1e-3,1e-3);
   
    cu = upperBounds(obj,2);
    cl = lowerBounds(obj,2);
    assertEqual(cu,0.002);
    assertEqual(cl,0);

function test_matRad_DoseConstraintFromObjective_computeConstrFunc
    objective = DoseObjectives.matRad_SquaredOverdosing(100,30);
    obj = DoseConstraints.matRad_DoseConstraintFromObjective(objective,1e-3,1e-3);

    dose = [25 28 30 23 35]';
    cDose = computeDoseConstraintFunction(obj,dose);
    cDose_expected = 5;
    assertElementsAlmostEqual(cDose,cDose_expected,'absolute',1e-4);

function test_matRad_DoseConstraintFromObjective_computeJacob
    objective = DoseObjectives.matRad_SquaredOverdosing(100,30);
    obj = DoseConstraints.matRad_DoseConstraintFromObjective(objective,1e-3,1e-3);

    dose = [25 28 30 23 35]';
    cDoseJacob  = computeDoseConstraintJacobian(obj,dose);
    cDoseJacob_expected = [0 0 0 0 2]';
    assertElementsAlmostEqual(cDoseJacob,cDoseJacob_expected,'absolute',1e-4);

function test_matRad_DoseConstraintFromObjective_getSet
    objective = DoseObjectives.matRad_SquaredOverdosing(100,30);
    obj = DoseConstraints.matRad_DoseConstraintFromObjective(objective,1e-3,1e-3);

    doseParams = getDoseParameters(obj);
    assertEqual(doseParams,30);
    
    obj = setDoseParameters(obj,60);
    assertEqual(obj.objective.parameters{1},60);

