function test_suite = test_directApertureOptimization
% The output should always be test_suite, and the function name the same as
% your file name

% To collect all tests defined below, this is needed in newer Matlab
% versions. test_functions will collect function handles to below test
% functions
test_functions = localfunctions();

% This will initialize the test suite, i.e., take the functions from
% test_functions, check if they contain "test", convert them into a MOxUnit
% Test Case, and add them to the test-runner
initTestSuite;

function [dij, cst, pln, sequencing] = helper_getSequencedData()
p = load('photons_testData.mat');
pln = p.pln;
pln.propSeq.sequencer   = 'siochi';
pln.propOpt.quantityOpt = 'physicalDose';
pln.propOpt.maxIter     = 1; % keep the test fast, we only check the interface
resultGUI  = matRad_sequencing(p.resultGUI, p.stf, pln);
dij        = p.dij;
cst        = p.cst;
sequencing = resultGUI.sequencing;

function test_dao_restoredApiAndFieldPreservation
if ~matRad_OptimizerIPOPT.isAvailable()
    moxunit_throw_test_skipped_exception('IPOPT optimizer not available!');
end
[dij, cst, pln, sequencing] = helper_getSequencedData();

% call with the (restored) 5-argument API, passing an optResult struct that
% already carries data which must be preserved
resultGUI.sequencing       = sequencing;   % has .beam and .apertureInfo
resultGUI.someExistingField = 42;

[resultGUI, optimizer] = matRad_directApertureOptimization(dij, cst, sequencing.apertureInfo, resultGUI, pln);

% dose cubes and weights are added
assertTrue(isfield(resultGUI, 'physicalDose'));
assertTrue(isfield(resultGUI, 'w'));
assertEqual(resultGUI.w, resultGUI.wDAO);

% aperture info stored at both the backward-compatible top-level location and
% under the sequencing field
assertTrue(isstruct(resultGUI.apertureInfo));
assertTrue(isstruct(resultGUI.sequencing.apertureInfo));

% pre-existing (also nested) data in the passed-in optResult is preserved
assertTrue(isfield(resultGUI.sequencing, 'beam'));
assertEqual(resultGUI.someExistingField, 42);

% the used optimizer object is returned
assertTrue(isa(optimizer, 'matRad_OptimizerIPOPT'));
