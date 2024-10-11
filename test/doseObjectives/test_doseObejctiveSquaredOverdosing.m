function test_suite = test_doseObejctiveSquaredOverdosing 
    
test_functions=localfunctions();

initTestSuite;

function test_doseObjective_SquaredOverdosing_construct
    obj = DoseObjectives.matRad_SquaredOverdosing();   %default
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseObjectives.matRad_SquaredOverdosing'));
    assertEqual(obj.parameters{1},30);
    assertEqual(obj.penalty,1);
      
    obj = DoseObjectives.matRad_SquaredOverdosing(100,60); 
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseObjectives.matRad_SquaredOverdosing'));

function test_doseObjective_SquaredOverdosing_computeObjFunction

    obj = DoseObjectives.matRad_SquaredOverdosing(100,60); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    expected_fDose = 5.8;
   
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);
    
    dose = [55,58,60,NaN,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    
    assertTrue(isnan(test_fDose));


function test_doseObjective_SquaredOverdosing_computeGradient

   obj = DoseObjectives.matRad_SquaredOverdosing(100,60); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveGradient(obj,dose);
    expected_fDose = [0         0         0    0.8000    2.0000]'; 
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);

    dose = [55,58,60,NaN,65]';
    test_fDose = computeDoseObjectiveGradient(obj,dose);
    
    assertTrue(isnan(test_fDose(4)));