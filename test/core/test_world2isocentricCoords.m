function test_suite = test_world2isocentricCoords

test_functions = localfunctions();

initTestSuite;


function test_empty_world2isocentricCoords
    ct = get_testCtHelper();
    assertExceptionThrown(@()matRad_world2isocentricCoords());
    assertExceptionThrown(@()matRad_world2isocentricCoords([],ct));


function test_basic_world2isocentricCoords
    ct = get_testCtHelper();
    expected = [3 36 45];
    assertEqual(matRad_world2isocentricCoords([-30 3 12], ct),expected);
%   out of bounds 
    assertExceptionThrown(@()matRad_world2isocentricCoords([-30,4,45],ct));


function ct = get_testCtHelper()
    ct.x = -30:3:27;
    ct.y = -30:3:27;
    ct.z = -30:3:27;
    
    ct.resolution.x = 3;
    ct.resolution.y = 3;
    ct.resolution.z = 3;
    
    ct.cubeDim = [20 20 20 ];