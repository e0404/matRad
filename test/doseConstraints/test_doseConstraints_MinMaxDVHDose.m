function test_suite = test_doseConstraints_MinMaxDVHDose

test_functions=localfunctions();

initTestSuite;

function test_DoseConstraints_matRad_MinMaxDVHDose
    obj = DoseConstraints.matRad_MinMaxDVH();   %default
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseConstraints.matRad_MinMaxDVH'));
    assertEqual(obj.parameters{1},30);
    assertEqual(obj.parameters{2},0);
    assertEqual(obj.parameters{3},100);
    
    obj = DoseConstraints.matRad_MinMaxDVH(10,30);
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseConstraints.matRad_MinMaxDVH'));

function test_DoseConstraints_matRad_DVHDose_Bounds
    obj = DoseConstraints.matRad_MinMaxDVH(30,0,100);
    cu = upperBounds(obj,2);
    cl = lowerBounds(obj,2);
    assertEqual(cu,1);
    assertEqual(cl,0);


function test_DoseConstraints_matRad_DVHDose_computeDoseConstr
    obj = DoseConstraints.matRad_MinMaxDVH(60,0,100);
    dose = [55 58 60 61 65]';
    cDose = computeDoseConstraintFunction(obj,dose);
    expected_cDose = 0.6;
    assertEqual(cDose,expected_cDose);

function test_DoseConstraints_matRad_DVHDose_computeJacobian
    obj = DoseConstraints.matRad_MinMaxDVH(60,0,100);
    dose = [55 58 60 61 65]';
    cDoseJacob  = computeDoseConstraintJacobian(obj,dose);
    expected_Jacob =[0.0018; 0.0218; 0.0460; 0.0375; 0.0018];
    assertElementsAlmostEqual(cDoseJacob,expected_Jacob,'absolute',1e-4);
    