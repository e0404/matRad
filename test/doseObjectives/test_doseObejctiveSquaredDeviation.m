function test_suite = test_doseObejctiveSquaredDeviation 
    
test_functions=localfunctions();

initTestSuite;

function test_doseObjective_squaredDeviation_construct
    
    obj = DoseObjectives.matRad_SquaredDeviation();   %default
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseObjectives.matRad_SquaredDeviation'));
    assertEqual(obj.parameters{1},60);
    assertEqual(obj.penalty,1);
      
    obj = DoseObjectives.matRad_SquaredDeviation(100,60); 
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseObjectives.matRad_SquaredDeviation'));

function test_doseObjective_SquaredDeviation_computeObjFunction

    obj = DoseObjectives.matRad_SquaredDeviation(100,60); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    expected_fDose = 11.6000;
   
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);
    
    dose = [55,58,60,NaN,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    
    assertTrue(isnan(test_fDose));


function test_doseObjective_SquaredDeviation_computeGradient

   obj = DoseObjectives.matRad_SquaredDeviation(100,60); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveGradient(obj,dose);
    expected_fDose = [-2.0000; -0.8000; 0; 0.8000; 2.0000]; 
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);

    dose = [55,58,60,NaN,65]';
    test_fDose = computeDoseObjectiveGradient(obj,dose);
    
    assertTrue(isnan(test_fDose(4)));