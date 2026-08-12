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

function test_daoRestoredApiAndFieldPreservation
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

function test_daoScaleToPrescription
if ~matRad_OptimizerIPOPT.isAvailable()
    moxunit_throw_test_skipped_exception('IPOPT optimizer not available!');
end
[dij, cst, pln, sequencing] = helper_getSequencedData();

% reference run without prescription scaling
resultRef = matRad_directApertureOptimization(dij, cst, sequencing.apertureInfo, struct(), pln);

% run with the plan scaled such that the target D95 reaches the prescription
pln.propOpt.scaleToPrescription  = true;
pln.propOpt.prescribedDose       = 45; % [Gy] over all fractions
pln.propOpt.prescriptionStructIx = find(strcmp(cst(:, 3), 'TARGET'));
resultGUI = matRad_directApertureOptimization(dij, cst, sequencing.apertureInfo, struct(), pln);

factor = resultGUI.apertureInfo.prescriptionScaleFactor;
assertTrue(isscalar(factor) && isfinite(factor) && factor > 0);

% weight vectors stay consistent (regression: wDao/wDAO field case mismatch)
assertEqual(resultGUI.w, resultGUI.wDAO);

% the whole plan is scaled consistently by the same factor
assertElementsAlmostEqual(resultGUI.w, resultRef.w * factor, 'relative', 1e-10);
assertElementsAlmostEqual(resultGUI.physicalDose, resultRef.physicalDose * factor, 'relative', 1e-10);

% the target D95 of the scaled plan matches the prescribed dose per fraction
qi = matRad_calcQualityIndicators(cst, pln, resultGUI.physicalDose);
assertElementsAlmostEqual(qi(pln.propOpt.prescriptionStructIx).D_95, ...
                          pln.propOpt.prescribedDose / pln.numOfFractions, 'relative', 1e-6);
