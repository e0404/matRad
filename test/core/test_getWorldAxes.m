function test_suite = test_getWorldAxes

test_functions = localfunctions();

initTestSuite;

function test_basic_getWorldAxes
    assertExceptionThrown(@()matRad_getWorldAxes());
   
    ct = get_testCtHelper();
    expectedX =  [-5 -4 -3 -2 -1 0 1 2 3 4];
    expectedY = [-10 -8 -6 -4 -2 0 2 4 6 8];
    expectedZ = [-15 -12 -9 -6 -3 0 3 6 9 12];
    resultCt =  matRad_getWorldAxes(ct);
    assertEqual(resultCt.x,expectedX);
    assertEqual(resultCt.y,expectedY);
    assertEqual(resultCt.z,expectedZ);
    
function ct = get_testCtHelper()
    ct.resolution.x = 1;
    ct.resolution.y = 2;
    ct.resolution.z = 3;
    ct.cubeDim = [10 10 10 ];