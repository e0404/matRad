## matRad Tests
matRad tests run with [MOxUnit](https://github.com/MOxUnit/MOxUnit). Coverage is generated with [MOcov](https://github.com/MOcov/MOcov).

matRad runs through (most of the) examples (each one counting as one test) on top of defined unitTests using MOxUnits Syntax.
Multiple Matlab versions as well as Octave are tested, all running on Ubuntu.

### Testing Locally
1. As MOxUnit and MOcov are integrated as *submodules* in the submodules/ folder of the repository, make sure your submodules are initialized and up to date. Many git tools do that automatically on cloning a repository, but you can do this manually by running `git submodule update --init` in the repository root. 
   
   If you didn't clone the repository and are working with the downloaded source, you can do the following:
   1. Clone MOxUnit into a *separate* folder (not within or in subfolders of the matRad root dir).
   2. Navigate into the MOxUnit code folder (within the MOxUnit) and call `moxunit_set_path`.
   3. Navigate back into the matRad root folder.
2. Run tests on matRads test folder by calling `matRad_runTests(folder,withCoverage)` from the root folder. The arguments `folder` and `withCoverage` are optional. If `folder` is not given, it will default to `test`. A specific set of tests can be run, then, for example, by calling `matRad_runTests('test/tools')` to run tests in the `test/tools` folder.

### Writing Tests
Here's how a basic Test File can look like, based on [test_version.m](tools/test_version.m)
```matlab
function test_suite = test_version
%The output should always be test_suite, and the function name the same as
%your file name

%% Header
% The header is required to be in this format for automatic test collection
% by MOxUnit

%To collect all tests defined below, this is needed in newer Matlab
%versions. test_functions will collect function handles to below test
%functions
test_functions=localfunctions(); 

% This will initialize the test suite, i.e., take the functions from
% test_functions, check if they contain "test", convert them into a MOxUnit
% Test Case, and add them to the test-runner
initTestSuite;

%% Custom Tests
% Tests use assert*-like Functions to check outputs etc:
% assertTrue(a) - a is true
% assertFalse(a) - a is false
% assertEqual(a,b) - a and be are equal (isequal)
% assertElementsAlmostEqual(a,b) - numerical test for all vector / matrix elements. Has Additional arguments for absolute / relative tolerance 
% assertVectorsAlmostEqual(a,b) - numerical test using vector norm
% assertExceptionThrown(f,id) - test if exception of id is thrown. Take care of Octave issues with exception id (or don't provide id)
% Check MOxUnit for more information or look at other tests

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
```

### Testing a specific engine, generator or model
If the dose engine, stf generator, sequencer, biological model or scenario model requested in `pln` cannot be used, matRad warns and falls back to the default.
With the reduced log level of the test runner, this warning is not shown, so a test meant to exercise a specific class may pass while it actually runs the default.
Such tests should enable strict selection, which turns the fallback into an error.
`helper_setStrictSelection` returns a cleanup object that restores the previous setting when the test ends:
```matlab
function test_myEngineCalculatesDose
cleanup = helper_setStrictSelection(true); %#ok<NASGU>
pln.propDoseCalc.engine = 'MyEngine';
...
```
To check warnings or other messages independent of the log level, `helper_captureLog` returns the log entries recorded while calling a function, and `helper_logContains` searches them:
```matlab
logEntries = helper_captureLog(@() DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln));
assertTrue(helper_logContains(logEntries, 'WARNING', 'not available'));
```
