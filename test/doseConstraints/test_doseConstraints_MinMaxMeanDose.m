function test_suite = test_doseConstraints_MinMaxMeanDose

test_functions=localfunctions();

initTestSuite;

function test_DoseConstraints_matRad_MinMaxMeanDose

    obj = DoseConstraints.matRad_MinMaxMeanDose();   %default
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseConstraints.matRad_MinMaxMeanDose'));
    assertEqual(obj.parameters{1},0);
    assertEqual(obj.parameters{2},30);
    
    obj = DoseConstraints.matRad_MinMaxMeanDose(10,30);
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseConstraints.matRad_MinMaxMeanDose'));


function test_DoseConstraints_matRad_MeanDose_Bounds
    obj = DoseConstraints.matRad_MinMaxMeanDose(10,30);
    cu = upperBounds(obj,2);
    cl = lowerBounds(obj,2);
    assertEqual(cu,30);
    assertEqual(cl,10);


function test_DoseConstraints_matRad_MeanDose_ConstrJacob
    obj = DoseConstraints.matRad_MinMaxMeanDose(10,30);
    dose = [55,58,60,62,65]';
    cDose = computeDoseConstraintFunction(obj,dose);
    assertEqual(cDose,60);

    cDoseJacob  = computeDoseConstraintJacobian(obj,dose);
    expected_cJaco = [0.2 0.2 0.2 0.2 0.2]';
    assertEqual(cDoseJacob,expected_cJaco);