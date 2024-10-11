function test_suite = test_doseConstraintsMinMaxDose

test_functions=localfunctions();

initTestSuite;

function test_DoseConstraints_matRad_MinMaxDose

    obj = DoseConstraints.matRad_MinMaxDose();   %default
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseConstraints.matRad_MinMaxDose'));
    assertEqual(obj.parameters{1},0);
    assertEqual(obj.parameters{2},30);
    assertEqual(obj.parameters{3},1);
    
    
    obj = DoseConstraints.matRad_MinMaxDose(10,30,'approx');
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseConstraints.matRad_MinMaxDose'));

function test_DoseConstraints_matRad_MinMaxDose_upperBounds
    obj = DoseConstraints.matRad_MinMaxDose(-5,Inf,'approx');
    assertTrue(isempty(upperBounds(obj,2)));
    
    obj = DoseConstraints.matRad_MinMaxDose(50,Inf,'approx');
    assertTrue(isinf(upperBounds(obj,2)));
    
    obj = DoseConstraints.matRad_MinMaxDose(-5,70,'approx');
    assertEqual(upperBounds(obj,2),70);
    
    obj = DoseConstraints.matRad_MinMaxDose(50,70,'approx');
    assertEqual(upperBounds(obj,2),[Inf , 70]');
    
    obj = DoseConstraints.matRad_MinMaxDose(50,70,'voxelwise');
    assertEqual(upperBounds(obj,2),[70 70]');


function test_DoseConstraints_matRad_MinMaxDose_lowerBounds
    obj = DoseConstraints.matRad_MinMaxDose(-5,Inf,'approx');
    assertTrue(isempty(lowerBounds(obj,2)));
    
    obj = DoseConstraints.matRad_MinMaxDose(50,Inf,'approx');
    assertEqual(lowerBounds(obj,2),50);
    
    obj = DoseConstraints.matRad_MinMaxDose(-5,70,'approx');
    assertEqual(lowerBounds(obj,2),0);
    
    obj = DoseConstraints.matRad_MinMaxDose(50,70,'approx');
    assertEqual(lowerBounds(obj,2),[50, 0]');
    
    obj = DoseConstraints.matRad_MinMaxDose(50,70,'voxelwise');
    assertEqual(lowerBounds(obj,2),[50 50]');


function test_DoseConstraints_MinMaxDose_getDoseConstrntJacbnStrct
    obj = DoseConstraints.matRad_MinMaxDose(-5,Inf,'approx');
    assertTrue(isempty(getDoseConstraintJacobianStructure(obj,2)));
    
    obj = DoseConstraints.matRad_MinMaxDose(50,70,'approx');
    assertEqual(size(getDoseConstraintJacobianStructure(obj,2)),[2,2]);

    obj = DoseConstraints.matRad_MinMaxDose(-5,70,'approx');
    assertEqual(size(getDoseConstraintJacobianStructure(obj,2)),[2,1]);
    
    obj = DoseConstraints.matRad_MinMaxDose(50,70,'voxelwise');
    assertEqual(size(getDoseConstraintJacobianStructure(obj,2)),[2,2]);



function test_DoseConstraints_MinMaxDose_computeConstrFunction
    
    obj = DoseConstraints.matRad_MinMaxDose(-5,Inf ,'approx');    
    dose = [55,58,60,62,65]';
    test_cDose = computeDoseConstraintFunction(obj,dose);
    assertTrue(isempty(test_cDose));

    obj = DoseConstraints.matRad_MinMaxDose(50,70,'approx');
    
    dose = [55,58,60,62,65]';
    test_cDose = computeDoseConstraintFunction(obj,dose);
    expected_cDose = [55; 65];
    
    assertElementsAlmostEqual(test_cDose,expected_cDose,'absolute', 1e-4);

    obj = DoseConstraints.matRad_MinMaxDose(50,Inf,'approx');
    
    dose = [55,58,60,62,65]';
    test_cDose = computeDoseConstraintFunction(obj,dose);
    expected_cDose = [55];

    assertElementsAlmostEqual(test_cDose,expected_cDose,'absolute', 1e-4);

    obj = DoseConstraints.matRad_MinMaxDose(-5 ,70,'approx');
    
    dose = [55,58,60,62,65]';
    test_cDose = computeDoseConstraintFunction(obj,dose);
    expected_cDose = [65];
    
    assertElementsAlmostEqual(test_cDose,expected_cDose,'absolute', 1e-4);
    
    
    obj = DoseConstraints.matRad_MinMaxDose(50,70,'voxelwise');
    
    dose = [55,58,60,62,65]';
    test_cDose = computeDoseConstraintFunction(obj,dose);
    expected_cDose = [55,58,60,62,65]';
    
    assertElementsAlmostEqual(test_cDose,expected_cDose,'absolute', 1e-4);
    

function test_DoseConstraints_MinMaxDose_computeDoseConstraintJacobian
    
    obj = DoseConstraints.matRad_MinMaxDose(-5,Inf ,'approx');    
    dose = [55,58,60,62,65]';
    test_Jaco = computeDoseConstraintJacobian(obj,dose);
    assertTrue(isempty(test_Jaco));


    obj = DoseConstraints.matRad_MinMaxDose(50,70,'approx');    
    dose = [55,58,60,62,65]';

    test_Jaco = computeDoseConstraintJacobian(obj,dose);
    expected_Jaco = [ 1     0;
                      0     0;
                      0     0;
                      0     0;
                      0     1 ];
    assertElementsAlmostEqual(test_Jaco,expected_Jaco,'absolute', 1e-4);

    obj = DoseConstraints.matRad_MinMaxDose(50,Inf,'approx');
    dose = [55,58,60,62,65]';

    test_Jaco = computeDoseConstraintJacobian(obj,dose);
    expected_Jaco = [ 1 0 0 0 0]';   
    assertElementsAlmostEqual(test_Jaco,expected_Jaco,'absolute', 1e-4);

    obj = DoseConstraints.matRad_MinMaxDose(-5 ,70,'approx');
    dose = [55,58,60,62,65]';

    test_Jaco = computeDoseConstraintJacobian(obj,dose);
    expected_Jaco = [ 0 0 0 0 1 ]';   
    assertElementsAlmostEqual(test_Jaco,expected_Jaco,'absolute', 1e-4);

    obj = DoseConstraints.matRad_MinMaxDose(50,70,'voxelwise');    
    dose = [55,58,60,62,65]';

    test_Jaco = computeDoseConstraintJacobian(obj,dose);
    expected_Jaco = speye(5);
    assertElementsAlmostEqual(test_Jaco,expected_Jaco,'absolute', 1e-4);


