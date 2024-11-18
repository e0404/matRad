function test_suite = test_doseObejctiveMaxDVH 
    
test_functions=localfunctions();

initTestSuite;

function test_doseObjective_MaxDVH_construct

    obj = DoseObjectives.matRad_MaxDVH();   %default
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseObjectives.matRad_MaxDVH'));
    assertEqual(obj.parameters{1},30);
    assertEqual(obj.parameters{2},95);
    assertEqual(obj.penalty,1);
      
    obj = DoseObjectives.matRad_MaxDVH(100,30,95); 
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseObjectives.matRad_MaxDVH'));

function test_doseObjective_MaxDVH_computeObjFunction

    obj = DoseObjectives.matRad_MaxDVH(100,60,50); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    expected_fDose = 0;
    
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);

    obj = DoseObjectives.matRad_MaxDVH(100,50,95); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    expected_fDose = 5;
   
   
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);
    
    dose = [55,58,60,NaN,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    
    assertTrue(isnan(test_fDose));


function test_doseObjective_MaxDVH_computeGradient

   obj = DoseObjectives.matRad_MaxDVH(100,50,95); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveGradient(obj,dose);
    expected_fDose = [ 2 0 0 0 0]'; 
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);

    dose = [55,58,60,NaN,65]';
    test_fDose = computeDoseObjectiveGradient(obj,dose);
    
    assertTrue(isnan(test_fDose(4)));