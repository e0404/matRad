function test_suite = test_doseObejctiveSquaredUnderdosing 
    
test_functions=localfunctions();

initTestSuite;

function test_doseObjective_SquaredUnderdosing_construct

    obj = DoseObjectives.matRad_SquaredUnderdosing();   %default
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseObjectives.matRad_SquaredUnderdosing'));
    assertEqual(obj.parameters{1},60);
    assertEqual(obj.penalty,1);
      
    obj = DoseObjectives.matRad_SquaredUnderdosing(100,60); 
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseObjectives.matRad_SquaredUnderdosing'));

function test_doseObjective_SquaredUnderdosing_computeObjFunction

    obj = DoseObjectives.matRad_SquaredUnderdosing(100,60); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    expected_fDose = 5.8;
   
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);
    
    dose = [55,58,60,nan,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    
    assertTrue(isnan(test_fDose));


function test_doseObjective_SquaredUnderdosing_computeGradient

   obj = DoseObjectives.matRad_SquaredUnderdosing(100,60); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveGradient(obj,dose);
    expected_fDose = [ -2.0000 -0.8000 0 0 0]'; 
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);

    dose = [55,58,60,nan,65]';
    test_fDose = computeDoseObjectiveGradient(obj,dose);
    
    assertTrue(isnan(test_fDose(4)));