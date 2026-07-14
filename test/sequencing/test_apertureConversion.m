function test_suite = test_apertureConversion
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

function [pln, stf, sequencing] = helper_getSequencedData()
p = load('photons_testData.mat');
pln = p.pln;
pln.propSeq.sequencer = 'siochi';
resultGUI = matRad_sequencing(p.resultGUI, p.stf, pln);
stf = p.stf;
sequencing = resultGUI.sequencing;

function test_aperture2collimation_setsCollimation
[pln, stf, sequencing] = helper_getSequencedData();

[plnOut, stfOut] = matRad_aperture2collimation(pln, stf, sequencing, sequencing.apertureInfo);

% pln is switched to field-based collimation
assertEqual(plnOut.propStf.bixelWidth, 'field');
assertTrue(isfield(plnOut.propStf, 'collimation'));
assertTrue(isscalar(plnOut.propStf.collimation.fieldWidth));
assertTrue(isscalar(plnOut.propStf.collimation.leafWidth));

% stf now carries shapes instead of beamlets, one field ray per beam
assertEqual(numel(stfOut), numel(stf));
for i = 1:numel(stfOut)
    assertEqual(stfOut(i).numOfRays, 1);
    assertTrue(isfield(stfOut(i).ray, 'shapes'));
    assertEqual(numel(stfOut(i).ray.shapes), sequencing.beam(i).numOfShapes);
end

function test_aperture2collimation_missingApertureInfoErrors
[pln, stf, sequencing] = helper_getSequencedData();
sequencing = rmfield(sequencing, 'apertureInfo');

% without apertureInfo in the struct and no 4th argument the call must error
assertExceptionThrown(@() matRad_aperture2collimation(pln, stf, sequencing), 'matRad:Error');
