function test_suite = test_world2cubeIndex

test_functions = localfunctions();

initTestSuite;


function test_empty_world2cubeIndex
    ct = helper_getTestCt();
    assertExceptionThrown(@()matRad_world2cubeIndex());
    assertExceptionThrown(@()matRad_world2cubeIndex([],ct));

function test_basic_world2cubeIndex
    ct = helper_getTestCt();
    wCoord = [-30 3 12];
    expected = [12 1 15];
    assertEqual(matRad_world2cubeIndex(wCoord, ct),expected);

    wCoord = [-30 3 12; 12 -30 3];
    expected = [12 1 15; 1 15 12];
    assertEqual(matRad_world2cubeIndex(wCoord, ct),expected);

%   out of bounds 
    assertExceptionThrown(@()matRad_world2cubeIndex([-30,4,35],ct));


function ct = helper_getTestCt()
    ct.x = -30:3:27;
    ct.y = -30:3:27;
    ct.z = -30:3:27;
    
    ct.resolution.x = 3;
    ct.resolution.y = 3;
    ct.resolution.z = 3;
    
    ct.cubeDim = [20 20 20 ];