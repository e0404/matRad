function test_suite = test_sampleDij

test_functions = localfunctions();

initTestSuite;

function [ix, dose, radDepth, radDistSq, sampling] = helper_samplingInput()
n = 2000;
ix = (101:100 + n)';
radDistSq = linspace(0, 60, n)'.^2;
radDepth = mod((0:n - 1)', 100) + 0.5;
dose = exp(-sqrt(radDistSq) / 10);
sampling = struct('relDoseThreshold', 0.01, 'latCutOff', 20, 'type', 'radius', 'deltaRadDepth', 5);

function test_sampleDijKeepsCore
[ix, dose, radDepth, radDistSq, sampling] = helper_samplingInput();
[ixNew, doseNew] = matRad_sampleDij(ix, dose, radDepth, radDistSq, sampling);
isCore = radDistSq < sampling.latCutOff^2;
assertEqual(ixNew(1:nnz(isCore)), ix(isCore));
assertEqual(doseNew(1:nnz(isCore)), dose(isCore));
assertTrue(numel(ixNew) < numel(ix));
assertEqual(numel(unique(ixNew)), numel(ixNew));
assertTrue(all(ismember(ixNew, ix)));

function test_sampleDijPreservesExpectedDose
[ix, dose, radDepth, radDistSq, sampling] = helper_samplingInput();
total = 0;
numRuns = 200;
for i = 1:numRuns
    [~, doseNew] = matRad_sampleDij(ix, dose, radDepth, radDistSq, sampling);
    total = total + sum(doseNew);
end
assertElementsAlmostEqual(total / numRuns, sum(dose), 'relative', 0.02);

function test_sampleDijBixelWidthWidensCore
[ix, dose, radDepth, radDistSq, sampling] = helper_samplingInput();
sampling.latCutOff = 55;
[ixNew, doseNew] = matRad_sampleDij(ix, dose, radDepth, radDistSq, sampling, 10);
assertEqual(ixNew, ix);
assertEqual(doseNew, dose);

function test_sampleDijDoseType
[ix, dose, radDepth, radDistSq, sampling] = helper_samplingInput();
sampling.type = 'dose';
[ixNew, doseNew] = matRad_sampleDij(ix, dose, radDepth, radDistSq, sampling);
isCore = dose > sampling.relDoseThreshold * max(dose);
assertEqual(doseNew(1:nnz(isCore)), dose(isCore));
assertTrue(numel(ixNew) <= numel(ix));

function test_sampleDijLogicalIndex
[ix, dose, radDepth, radDistSq, sampling] = helper_samplingInput();
mask = false(max(ix), 1);
mask(ix) = true;
sampling.latCutOff = 100;
ixNew = matRad_sampleDij(mask, dose, radDepth, radDistSq, sampling);
assertEqual(ixNew, ix);

function test_sampleDijUnknownType
[ix, dose, radDepth, radDistSq, sampling] = helper_samplingInput();
sampling.type = 'unknown';
assertExceptionThrown(@() matRad_sampleDij(ix, dose, radDepth, radDistSq, sampling));
