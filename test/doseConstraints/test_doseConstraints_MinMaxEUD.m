function test_suite = test_doseConstraints_MinMaxEUD

test_functions=localfunctions();

initTestSuite;

function test_DoseConstraints_matRad_MinMaxEUD
    obj = DoseConstraints.matRad_MinMaxEUD();   %default
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseConstraints.matRad_MinMaxEUD'));
    assertEqual(obj.parameters{1},5);
    assertEqual(obj.parameters{2},0);
    assertEqual(obj.parameters{3},30);
    
    obj = DoseConstraints.matRad_MinMaxEUD(2,0,30);
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseConstraints.matRad_MinMaxEUD'));

function test_DoseConstraints_matRad_EUDDose_Bounds
    obj = DoseConstraints.matRad_MinMaxEUD(5,0,30);
    cu = upperBounds(obj,2);
    cl = lowerBounds(obj,2);
    assertEqual(cu,30);
    assertEqual(cl,0);

function test_DoseConstraints_matRad_EUDDose_computeDoseConstr
    obj = DoseConstraints.matRad_MinMaxEUD(5,0,30);
    dose = [55 58 60 62 65]';
    cDose = computeDoseConstraintFunction(obj,dose);
    expected_cDose = 60.3829;
    assertElementsAlmostEqual(cDose,expected_cDose,'absolute',1e-4);

function test_DoseConstraints_matRad_EUDDose_computeJacobian
    obj = DoseConstraints.matRad_MinMaxEUD(5,0,60);
    dose = [55 58 60 61 65]';
    cDoseJacob  = computeDoseConstraintJacobian(obj,dose);
    expected_Jacob =[0.1397; 0.1727; 0.1978; 0.2113; 0.2724];
    assertElementsAlmostEqual(cDoseJacob,expected_Jacob,'absolute',1e-4);
