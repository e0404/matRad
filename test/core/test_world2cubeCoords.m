function test_suite = test_world2cubeCoords

test_functions = localfunctions();

initTestSuite;


function test_empty_world2cubeCoords
    ct = get_testCtHelper();
    assertExceptionThrown(@()matRad_world2cubeCoords());
    assertExceptionThrown(@()matRad_world2cubeCoords([],ct));


function test_basic_world2cubeCoords
    ct = get_testCtHelper();
    expected = [1 12 15];
    assertEqual(matRad_world2cubeCoords([-30 3 12], ct),expected);
%   out of bounds 
    assertExceptionThrown(@()matRad_world2cubeCoords([-30,4,35],ct));


function ct = get_testCtHelper()
    ct.x = -30:3:27;
    ct.y = -30:3:27;
    ct.z = -30:3:27;
    
    ct.resolution.x = 3;
    ct.resolution.y = 3;
    ct.resolution.z = 3;
    
    ct.cubeDim = [20 20 20 ];