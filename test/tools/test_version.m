function test_suite = test_version

test_functions=localfunctions();

initTestSuite;

function test_matrad_version
    assertTrue(ischar(matRad_version()));
    [str,full] = matRad_version();
    assertTrue(isstruct(full));
    assertTrue(ischar(str));

function test_matrad_environment
    env = matRad_getEnvironment();
    assertTrue(ischar(env));
    [~,verstr] = matRad_getEnvironment();
    assertTrue(ischar(verstr));
    


