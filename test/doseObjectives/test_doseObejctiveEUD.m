function test_suite = test_doseObejctiveEUD 
    
test_functions=localfunctions();

initTestSuite;

function test_doseObjective_EUD_construct
    
    obj = DoseObjectives.matRad_EUD();   %default
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseObjectives.matRad_EUD'));
    assertEqual(obj.parameters{1},0);
    assertEqual(obj.parameters{2},3.5);
    assertEqual(obj.penalty,1);
      
    obj = DoseObjectives.matRad_EUD(100,60,3.5); 
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseObjectives.matRad_EUD'));

function test_doseObjective_EUD_computeObjFunction

    obj = DoseObjectives.matRad_EUD(100,60,3.5); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    expected_fDose = 0.0579;
   
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);
    
    dose = [55,58,60,NaN,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    
    assertTrue(isnan(test_fDose));


function test_doseObjective_EUD_computeGradient

   obj = DoseObjectives.matRad_EUD(100,10,2); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveGradient(obj,dose);
    expected_fDose = [ 18.3392   19.3395   20.0064   20.6733   21.6736]'; 
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);

    dose = [55,58,60,NaN,65]';
    assertExceptionThrown(@()computeDoseObjectiveGradient(obj,dose),'matRad:Error')
  