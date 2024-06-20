function test_suite = test_runBrachyTestSuite

test_functions=localfunctions();

initTestSuite;


function test_setup()
    if ~exist('matRad_cfg', 'var')
        addpath(fileparts(fileparts(pwd)));
        matRad_cfg = MatRad_Config.instance();
        matRad_rc;
    end


function test_stf()
    test_setup();  % Ensure setup is called
    testStf = test_generateBrachyStfTest;
    testResultStf = run(testStf);
    assert(~isempty(testResultStf), 'test_stf: Test results should not be empty');


function test_calcDose()
    test_setup();  % Ensure setup is called
    testCalcDose = test_calcBrachyDoseTest;
    testResultCalcDose = run(testCalcDose);
    assert(~isempty(testResultCalcDose), 'test_calcDose: Test results should not be empty');


