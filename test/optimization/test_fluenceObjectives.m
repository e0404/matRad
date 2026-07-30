function test_suite = test_fluenceObjectives
% Unit tests for the fluence domain optimization functions, i.e. objectives
% acting on the bixel weight vector directly instead of on a dose related
% quantity, and for the helpers collecting them from a plan.
%
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

function dij = helper_syntheticDij(beamSizes)
% Minimal dij carrying only what a fluence function needs: the mapping from
% bixel to beam. No geometry, no dose influence.
dij.numOfBeams       = numel(beamSizes);
dij.numOfRaysPerBeam = beamSizes;
dij.totalNumOfBixels = sum(beamSizes);

dij.beamNum = zeros(sum(beamSizes), 1);
offset = 0;
for b = 1:numel(beamSizes)
    dij.beamNum(offset + (1:beamSizes(b))) = b;
    offset = offset + beamSizes(b);
end

function w = helper_testFluence(dij, seed)
% a reproducible modulated fluence
if nargin < 2
    seed = 11;
end
rng(seed);
w = 1 + 0.7 * rand(dij.totalNumOfBixels, 1);

function f = helper_referenceRelative(w, dij, beams)
% mean over beams of (within-beam variance / within-beam mean square)
f = 0;
for b = beams
    wb = w(dij.beamNum == b);
    f = f + (mean((wb - mean(wb)).^2) / mean(wb.^2)) / numel(beams);
end

function f = helper_referenceAbsolute(w, dij, beams)
% pooled within-beam variance, normalized by the number of active bixels
numActive = nnz(ismember(dij.beamNum, beams));
f = 0;
for b = beams
    wb = w(dij.beamNum == b);
    f = f + sum((wb - mean(wb)).^2) / numActive;
end

function test_fluenceObjectiveClassHierarchy
% A fluence objective must not be a dose optimization function - no dose is
% involved - but both must share the common ancestor that carries the
% parameter/serialization interface.
fObj = FluenceObjectives.matRad_FluenceVariance(3);

assertTrue(isa(fObj, 'matRad_OptimizationFunction'));
assertTrue(isa(fObj, 'matRad_FluenceOptimizationFunction'));
assertFalse(isa(fObj, 'matRad_DoseOptimizationFunction'));
assertEqual(fObj.penalty, 3);
assertEqual(fObj.getNormalization(), 'relative');

% reparenting must not have disturbed the dose objectives
doseObj = DoseObjectives.matRad_SquaredDeviation(100, 60);
assertTrue(isa(doseObj, 'matRad_OptimizationFunction'));
assertTrue(isa(doseObj, 'matRad_DoseOptimizationFunction'));
assertEqual(doseObj.robustness, 'none');
assertElementsAlmostEqual(doseObj.getDoseParameters(), 60);

% including the legacy struct format
legacy = matRad_DoseOptimizationFunction.createInstanceFromStruct( ...
                                                                  struct('type', 'square deviation', 'penalty', 5, 'dose', 42));
assertTrue(isa(legacy, 'DoseObjectives.matRad_SquaredDeviation'));
assertElementsAlmostEqual(legacy.getDoseParameters(), 42);

function test_fluenceObjectiveSetupAndCompatibility
% The objective is a pure specification until it is set up for a dij.
dij = helper_syntheticDij([35, 20, 45]);
fObj = FluenceObjectives.matRad_FluenceVariance(1);

[isCompat, msg] = fObj.isCompatible(dij.totalNumOfBixels);
assertFalse(isCompat);
assertFalse(isempty(msg));

fObj = fObj.setupForDij(dij);
assertEqual(fObj.numOfBixels, dij.totalNumOfBixels);
assertEqual(numel(fObj.bixelIx), dij.totalNumOfBixels);
assertTrue(fObj.isCompatible(dij.totalNumOfBixels));
assertFalse(fObj.isCompatible(dij.totalNumOfBixels + 1));

% evaluating without a setup must fail rather than return nonsense
fresh = FluenceObjectives.matRad_FluenceVariance(1);
assertExceptionThrown(@() fresh.computeFluenceObjectiveFunction(ones(dij.totalNumOfBixels, 1)));

function test_fluenceVarianceAbsoluteIsTheQuadraticForm
% 'absolute' must be exactly w'Sw for the S the static builder returns, and S
% must be symmetric, positive semidefinite and block diagonal by beam.
beamSizes = [35, 20, 45];
dij = helper_syntheticDij(beamSizes);
w = helper_testFluence(dij);

fObj = FluenceObjectives.matRad_FluenceVariance(1, 'absolute');
fObj = fObj.setupForDij(dij);

fVal = fObj.computeFluenceObjectiveFunction(w);
assertElementsAlmostEqual(fVal, helper_referenceAbsolute(w, dij, 1:3), 'absolute', 1e-12);

quadForm = FluenceObjectives.matRad_FluenceVariance.getQuadraticForm(dij);
assertElementsAlmostEqual(w' * quadForm * w, fVal, 'absolute', 1e-12);

quadFull = full(quadForm);
assertElementsAlmostEqual(quadFull, quadFull', 'absolute', 1e-12);
assertTrue(min(eig(quadFull)) > -1e-10, 'quadratic form is not positive semidefinite');

% no coupling between beams: a flat fluence per beam is free even when the
% beams sit at completely different levels
assertEqual(nnz(quadForm(1:beamSizes(1), beamSizes(1) + 1:end)), 0);
wFlat = [ones(beamSizes(1), 1); 7 * ones(beamSizes(2), 1); 3 * ones(beamSizes(3), 1)];
assertElementsAlmostEqual(fObj.computeFluenceObjectiveFunction(wFlat), 0, 'absolute', 1e-12);

function test_fluenceVarianceRelativeIsNormalizedPerBeam
% 'relative' is the mean over beams of the squared coefficient of variation.
dij = helper_syntheticDij([35, 20, 45]);
w = helper_testFluence(dij);

fObj = FluenceObjectives.matRad_FluenceVariance(1, 'relative');
fObj = fObj.setupForDij(dij);

assertElementsAlmostEqual(fObj.computeFluenceObjectiveFunction(w), ...
                          helper_referenceRelative(w, dij, 1:3), 'absolute', 1e-12);

wFlat = [ones(35, 1); 7 * ones(20, 1); 3 * ones(45, 1)];
assertElementsAlmostEqual(fObj.computeFluenceObjectiveFunction(wFlat), 0, 'absolute', 1e-12);

function test_fluenceVarianceRelativeIgnoresBeamLevels
% Regression: the denominator used to be the mean square pooled over all
% beams. Since
%     sum_i w_i^2 = sum_b sum_{i in b} (w_i - m_b)^2 + sum_b n_b m_b^2 ,
% that put the differences between the beam levels into the denominator, so
% the objective could be lowered by driving the beams apart in intensity
% instead of by smoothing - on a VMAT arc it concentrated the plan on a few
% control points. Normalizing each beam with its own mean square makes the
% objective invariant under rescaling individual beams.
dij = helper_syntheticDij([35, 20, 45]);
w = helper_testFluence(dij);

fObj = FluenceObjectives.matRad_FluenceVariance(1, 'relative');
fObj = fObj.setupForDij(dij);
fBase = fObj.computeFluenceObjectiveFunction(w);

wSpread = w;
wSpread(dij.beamNum == 1) = 3.0 * wSpread(dij.beamNum == 1);
wSpread(dij.beamNum == 2) = 0.2 * wSpread(dij.beamNum == 2);
assertElementsAlmostEqual(fObj.computeFluenceObjectiveFunction(wSpread), fBase, 'relative', 1e-10, ...
                          'rescaling individual beams changed the objective');

% and it stays invariant under a global rescaling, so penalties transfer
assertElementsAlmostEqual(fObj.computeFluenceObjectiveFunction(10 * w), fBase, 'relative', 1e-10);

% the absolute variant is not scale invariant, by definition
fAbs = FluenceObjectives.matRad_FluenceVariance(1, 'absolute');
fAbs = fAbs.setupForDij(dij);
assertElementsAlmostEqual(fAbs.computeFluenceObjectiveFunction(10 * w), ...
                          100 * fAbs.computeFluenceObjectiveFunction(w), 'relative', 1e-10);

function test_fluenceVarianceGradients
% finite difference check of the analytic gradients
dij = helper_syntheticDij([12, 8, 15]);
w = helper_testFluence(dij, 7);

for normalization = {'relative', 'absolute'}
    fObj = FluenceObjectives.matRad_FluenceVariance(1, normalization{1});
    fObj = fObj.setupForDij(dij);

    gradAnalytic = fObj.computeFluenceObjectiveGradient(w);
    assertEqual(size(gradAnalytic), [dij.totalNumOfBixels, 1]);

    gradNumeric = zeros(size(w));
    h = 1e-6;
    for k = 1:numel(w)
        wPlus = w;
        wPlus(k) = wPlus(k) + h;
        wMinus = w;
        wMinus(k) = wMinus(k) - h;
        gradNumeric(k) = (fObj.computeFluenceObjectiveFunction(wPlus) - ...
                          fObj.computeFluenceObjectiveFunction(wMinus)) / (2 * h);
    end

    relErr = norm(gradAnalytic - gradNumeric) / norm(gradNumeric);
    assertTrue(relErr < 1e-6, sprintf('%s gradient off by %.3e', normalization{1}, relErr));
end

function test_fluenceVarianceHandlesBeamWithoutFluence
% A beam that carries no fluence has no modulation to speak of and must not
% produce NaN or Inf in either the value or the gradient.
dij = helper_syntheticDij([20, 15, 20]);
w = helper_testFluence(dij);
w(dij.beamNum == 2) = 0;

fObj = FluenceObjectives.matRad_FluenceVariance(1, 'relative');
fObj = fObj.setupForDij(dij);

assertTrue(isfinite(fObj.computeFluenceObjectiveFunction(w)));
assertTrue(all(isfinite(fObj.computeFluenceObjectiveGradient(w))));

assertTrue(isfinite(fObj.computeFluenceObjectiveFunction(zeros(size(w)))));
assertTrue(all(isfinite(fObj.computeFluenceObjectiveGradient(zeros(size(w))))));

function test_fluenceVarianceBeamSelection
% By default a VMAT plan is smoothed only where fluence is optimized, i.e. on
% the FMO beams; the remaining control points are bounded to zero anyway.
beamSizes = [35, 20, 45];
dij = helper_syntheticDij(beamSizes);
dijVMAT = dij;
dijVMAT.isFMOBeam = [true, false, true];

fObj = FluenceObjectives.matRad_FluenceVariance(1);
fObj = fObj.setupForDij(dijVMAT);

assertEqual(numel(fObj.bixelIx), beamSizes(1) + beamSizes(3));
inNonFMOBeam = fObj.bixelIx > beamSizes(1) & fObj.bixelIx <= beamSizes(1) + beamSizes(2);
assertFalse(any(inNonFMOBeam));

% an explicit beam selection overrides the default
fBeam = FluenceObjectives.matRad_FluenceVariance(1);
fBeam.beams = 2;
fBeam = fBeam.setupForDij(dij);
assertEqual(numel(fBeam.bixelIx), beamSizes(2));

% several bixels per ray (particles) need no special handling
dijIons = helper_syntheticDij([30, 20]);
dijIons.numOfRaysPerBeam = [10, 10];
fIons = FluenceObjectives.matRad_FluenceVariance(1);
fIons = fIons.setupForDij(dijIons);
assertEqual(numel(fIons.bixelIx), 50);

function test_fluenceObjectiveSerialization
% The objective must survive a struct round trip as a specification; the
% beamlet partition is rebuilt from the dij afterwards.
dij = helper_syntheticDij([35, 20, 45]);
w = helper_testFluence(dij);

fObj = FluenceObjectives.matRad_FluenceVariance(4, 'absolute');
fObj = fObj.setupForDij(dij);

s = struct(fObj);
assertEqual(s.className, 'FluenceObjectives.matRad_FluenceVariance');
assertEqual(s.penalty, 4);

restored = matRad_OptimizationFunction.createInstanceFromStruct(s);
assertTrue(isa(restored, 'FluenceObjectives.matRad_FluenceVariance'));
assertEqual(restored.penalty, 4);
assertEqual(restored.getNormalization(), 'absolute');
assertFalse(restored.isCompatible(dij.totalNumOfBixels));

restored = restored.setupForDij(dij);
assertElementsAlmostEqual(restored.computeFluenceObjectiveFunction(w), ...
                          fObj.computeFluenceObjectiveFunction(w), 'absolute', 1e-12);

function test_fluenceObjectiveRejectsInvalidNormalization
assertExceptionThrown(@() FluenceObjectives.matRad_FluenceVariance(1, 'quadratic'));

function test_getFluenceObjectivesCollectsFromPlnAndCst
% Fluence objectives are accepted both in pln.propOpt and next to the dose
% objectives of a structure in the cst.
dij = helper_syntheticDij([35, 20, 45]);

cst = cell(2, 6);
cst{1, 3} = 'TARGET';
cst{1, 6} = {DoseObjectives.matRad_SquaredDeviation(100, 60)};
cst{2, 3} = 'OAR';
cst{2, 6} = {DoseObjectives.matRad_SquaredOverdosing(50, 30), ...
             FluenceObjectives.matRad_FluenceVariance(7)};

pln.propOpt.fluenceObjectives = {FluenceObjectives.matRad_FluenceVariance(3)};

collected = matRad_getFluenceObjectives(pln, cst, dij);
assertEqual(numel(collected), 2);
assertElementsAlmostEqual(sort([collected{1}.penalty, collected{2}.penalty]), [3 7]);
for i = 1:numel(collected)
    assertTrue(collected{i}.isCompatible(dij.totalNumOfBixels));
end

% a bare object instead of a cell array is accepted too
plnBare.propOpt.fluenceObjectives = FluenceObjectives.matRad_FluenceVariance(1);
assertEqual(numel(matRad_getFluenceObjectives(plnBare, cell(0, 6), dij)), 1);

% passing no dij collects the specifications without setting them up, which
% is what direct aperture optimization uses to warn about them
unsetUp = matRad_getFluenceObjectives(pln, cst, []);
assertEqual(numel(unsetUp), 2);
assertFalse(unsetUp{1}.isCompatible(dij.totalNumOfBixels));

% plans without any fluence objective come back empty
assertTrue(isempty(matRad_getFluenceObjectives(struct(), cst(1, :), dij)));

function test_fractionateCstFunctionsSkipsFluenceObjectives
% The per-fraction rescaling applies to dose parameters only - a fluence
% objective has none and must come through untouched.
cst = cell(1, 6);
cst{1, 3} = 'TARGET';
cst{1, 6} = {DoseObjectives.matRad_SquaredDeviation(100, 60), ...
             FluenceObjectives.matRad_FluenceVariance(7)};

cst = matRad_fractionateCstFunctions(cst, 30);

assertElementsAlmostEqual(cst{1, 6}{1}.getDoseParameters(), 2);
assertTrue(isa(cst{1, 6}{2}, 'FluenceObjectives.matRad_FluenceVariance'));
assertEqual(cst{1, 6}{2}.penalty, 7);

% the legacy struct array format is still accepted
cstLegacy = cell(1, 6);
cstLegacy{1, 3} = 'TARGET';
cstLegacy{1, 6} = struct('type', 'square deviation', 'penalty', 800, 'dose', 60);
cstLegacy = matRad_fractionateCstFunctions(cstLegacy, 30);
assertTrue(isa(cstLegacy{1, 6}{1}, 'DoseObjectives.matRad_SquaredDeviation'));
assertElementsAlmostEqual(cstLegacy{1, 6}{1}.getDoseParameters(), 2);
