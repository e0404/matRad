function test_suite = test_doseObejctiveMeanDose 
    
test_functions=localfunctions();

initTestSuite;

function test_doseObjective_squaredDeviation_construct
    
    obj = DoseObjectives.matRad_MeanDose();   %default
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseObjectives.matRad_MeanDose'));
    assertEqual(obj.parameters{1},0);
    assertEqual(obj.parameters{2},1);
    assertEqual(obj.penalty,1);
      
    obj = DoseObjectives.matRad_MeanDose(100,0,'Linear'); 
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseObjectives.matRad_MeanDose'));

    obj = DoseObjectives.matRad_MeanDose(100,10,4);
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseObjectives.matRad_MeanDose'));
    assertEqual(obj.parameters{1},10);
    assertEqual(obj.parameters{2},1);

function test_doseObjective_MeanDose_computeObjFunction

    obj = DoseObjectives.matRad_MeanDose(100,10,'Linear'); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    expected_fDose = 50;
   
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);

    obj = DoseObjectives.matRad_MeanDose(100,0,'Quadratic'); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    expected_fDose = 3600;
   
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);
    
    dose = [55,58,60,NaN,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    
    assertTrue(isnan(test_fDose));


function test_doseObjective_MeanDose_computeGradient

   obj = DoseObjectives.matRad_MeanDose(100,0,'Linear'); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveGradient(obj,dose);
    expected_fDose = [0.2; 0.2; 0.2; 0.2; 0.2]; 
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);

    obj = DoseObjectives.matRad_MeanDose(100,0,'Quadratic'); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveGradient(obj,dose);
    expected_fDose = [24; 24; 24; 24; 24]; 
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);

    dose = [55,58,60,NaN,65]';
    test_fDose = computeDoseObjectiveGradient(obj,dose);
    
    assertTrue(isnan(test_fDose(4)));
