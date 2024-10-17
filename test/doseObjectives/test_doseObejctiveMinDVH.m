function test_suite = test_doseObejctiveMinDVH 
    
test_functions=localfunctions();

initTestSuite;

function test_doseObjective_MinDVH_construct

    obj = DoseObjectives.matRad_MinDVH();   %default
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseObjectives.matRad_MinDVH'));
    assertEqual(obj.parameters{1},60);
    assertEqual(obj.parameters{2},95);
    assertEqual(obj.penalty,1);
      
    obj = DoseObjectives.matRad_MinDVH(100,60,95); 
    assertTrue(isobject(obj));
    assertTrue(isa(obj, 'DoseObjectives.matRad_MinDVH'));

function test_doseObjective_MinDVH_computeObjFunction

    obj = DoseObjectives.matRad_MinDVH(100,60,50); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    expected_fDose = 0;
    
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);

    obj = DoseObjectives.matRad_MinDVH(100,60,95); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    expected_fDose = 5.8;
   
   
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);
    
    dose = [55,58,60,NaN,65]';
    test_fDose = computeDoseObjectiveFunction(obj,dose);
    
    assertTrue(isnan(test_fDose));


function test_doseObjective_MinDVH_computeGradient

   obj = DoseObjectives.matRad_MinDVH(100,60); 

    dose = [55,58,60,62,65]';
    test_fDose = computeDoseObjectiveGradient(obj,dose);
    expected_fDose = [ -2.0000 -0.8000 0 0 0]'; 
    assertElementsAlmostEqual(test_fDose,expected_fDose,'absolute', 1e-4);

    dose = [55,58,60,NaN,65]';
    test_fDose = computeDoseObjectiveGradient(obj,dose);
    
    assertTrue(isnan(test_fDose(4)));