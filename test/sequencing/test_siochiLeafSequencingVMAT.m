function test_suite = test_siochiLeafSequencingVMAT
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

function [stf, pln] = helper_getVmatStf(continuousAperture, gantryAngles, couchAngles, ...
                                        maxGantryAngleSpacing, maxDAOGantryAngleSpacing, maxFMOGantryAngleSpacing)
if nargin < 1
    continuousAperture = false;
end
if nargin < 2
    gantryAngles = [-180, 180];
end
if nargin < 3
    couchAngles = [0, 0];
end
if nargin < 4
    maxGantryAngleSpacing = 15;
end
if nargin < 5
    maxDAOGantryAngleSpacing = 30;
end
if nargin < 6
    maxFMOGantryAngleSpacing = 45;
end
p = load('photons_testData.mat', 'ct', 'cst', 'pln');
pln = p.pln;

pln.propStf.generator                = 'PhotonVMAT';
pln.propStf.gantryAngles             = gantryAngles;
pln.propStf.couchAngles              = couchAngles;
pln.propStf.maxGantryAngleSpacing    = maxGantryAngleSpacing;
pln.propStf.maxDAOGantryAngleSpacing = maxDAOGantryAngleSpacing;
pln.propStf.maxFMOGantryAngleSpacing = maxFMOGantryAngleSpacing;
pln.propStf.isoCenter                = matRad_getIsoCenter(p.cst, p.ct, 0);

pln.propSeq.continuousAperture = continuousAperture;
pln.propOpt.runVMAT = true;

stf = matRad_generateStf(p.ct, p.cst, pln);

function test_vmatAcceptsScalarCouchAngle
[stf, ~] = helper_getVmatStf(false, [-180, 180], 17);

assertTrue(~isempty(stf));
assertEqual([stf.couchAngle], 17 * ones(1, numel(stf)));

function test_vmatSupportsReverseDirectionArc
[stf, pln] = helper_getVmatStf(true, [180, -180], 0);

assertTrue(~isempty(stf));
assertTrue(all(diff([stf.gantryAngle]) < 0));
doseWidths = arrayfun(@(ix) stf(ix).arc.doseAngleBordersDiff, 1:numel(stf));
assertTrue(all(doseWidths > 0));
daoBeams = arrayfun(@(ix) stf(ix).arc.isDAOBeam, 1:numel(stf));
daoWidths = arrayfun(@(ix) stf(ix).arc.DAOAngleBordersDiff, find(daoBeams));
assertTrue(all(daoWidths > 0));

sequencer = helper_getSequencer(pln);
w = ones(sum([stf.numOfRays]), 1);
apertureInfo = sequencer.sequence(w, stf).apertureInfo;
assertTrue(all(isfinite(apertureInfo.bixelWeights)));
assertTrue(all(apertureInfo.bixelWeights >= 0));

p = load('photons_testData.mat', 'ct', 'cst');
[stfFine, apertureInfoFine] = matRad_refineApertureArc(p.ct, p.cst, pln, apertureInfo, 7.5);
assertTrue(all(diff([stfFine.gantryAngle]) < 0));
assertTrue(all(isfinite(apertureInfoFine.bixelWeights)));
assertTrue(any(apertureInfoFine.bixelWeights > 0));

function helper_assertTimeFracPartition(stf)
% The DAO angle border splits an interpolated beam's dose sector into a part
% whose time is taken from the last and one taken from the next DAO beam, so
% the two fractions must sum to exactly 1 -- matRad_daoVec2ApertureInfo feeds
% them into 1 ./ (fracLast ./ rotLast + fracNext ./ rotNext).
interpIx = find(arrayfun(@(s) ~s.arc.isFMOBeam && ~s.arc.isDAOBeam, stf));
assertTrue(~isempty(interpIx));
for ix = interpIx
    fracSum = stf(ix).arc.timeFracFromLastDAO + stf(ix).arc.timeFracFromNextDAO;
    assertElementsAlmostEqual(fracSum, 1, 'absolute', 1e-12, ...
                              sprintf('beam %d: time fractions sum to %g', ix, fracSum));
end

function test_vmatTimeFractionsPartitionDoseSector
% Regression: taking abs() of the numerators defeated the [0, 1] clamp for a
% DAO border lying outside the dose sector -- both fractions then clamped to
% 1 instead of to 0 and 1. Needs >= 3 dose angles per DAO sector (7.5 vs 30
% deg) to expose it; the default 15/30 fixture puts exactly one interpolated
% beam per sector, whose sector always contains the border.
[stfFwd, ~] = helper_getVmatStf(true, [-180, 180], 0, 7.5);
helper_assertTimeFracPartition(stfFwd);

[stfRev, ~] = helper_getVmatStf(true, [180, -180], 0, 7.5);
helper_assertTimeFracPartition(stfRev);

% both clamp ends must still be reachable, i.e. the test data really does
% contain beams whose sector does not straddle the DAO border
allFracs = arrayfun(@(s) s.arc.timeFracFromLastDAO, ...
                    stfFwd(arrayfun(@(s) ~s.arc.isFMOBeam && ~s.arc.isDAOBeam, stfFwd)));
assertTrue(any(allFracs == 0) || any(allFracs == 1));

function test_vmatSequencesArcsWithInterpolatedBeams
% 4/10/30 puts 3 dose angles in every DAO sector, so most beams carry no
% aperture of their own and the arc opens and closes with interpolated
% beams. Regression: leafTouching probed apertureInfo.beam(1).shape(1) to
% decide its sampling mode, which only exists on DAO beams - beam 1 is no
% longer guaranteed to be one. Configurations whose dose angles coincide
% with the DAO control points (e.g. 15/30/45) never exercise this.
for continuousAperture = [false true]
    [stf, pln] = helper_getVmatStf(continuousAperture, [-180, 180], [0, 0], 4, 10, 30);

    isDAO = arrayfun(@(s) s.arc.isDAOBeam, stf);
    assertTrue(~isDAO(1) && ~isDAO(end), 'fixture must open and close with interpolated beams');
    assertTrue(nnz(~isDAO) > nnz(isDAO), 'fixture must be dominated by interpolated beams');

    sequencer = helper_getSequencer(pln);
    w = ones(sum([stf.numOfRays]), 1);
    apertureInfo = sequencer.sequence(w, stf).apertureInfo;

    assertEqual(apertureInfo.totalNumOfShapes, nnz(isDAO));
    assertTrue(all(isfinite(apertureInfo.bixelWeights)));
    assertTrue(all(apertureInfo.bixelWeights >= 0));
    assertTrue(any(apertureInfo.bixelWeights > 0));

    % every beam - interpolated ones included - ends up with usable leaf
    % positions once matRad_daoVec2ApertureInfo has run
    for i = 1:numel(stf)
        shape = apertureInfo.beam(i).shape(1);
        assertTrue(all(isfinite(shape.leftLeafPos)) && all(isfinite(shape.rightLeafPos)));
        assertTrue(all(shape.leftLeafPos <= shape.rightLeafPos));
    end
end

function sequencer = helper_getSequencer(pln)
sequencer = matRad_SequencingPhotonsSiochiLeaf(pln);
sequencer.runVMAT = true;
sequencer.numLevels = 5;
sequencer.weightToMU = 100;

function test_siochiVMATIsVMATCapable
[~, pln] = helper_getVmatStf();
sequencer = helper_getSequencer(pln);
assertTrue(isa(sequencer, 'matRad_PhotonSequencerVMATAbstract'));
assertTrue(sequencer.isVMATCapable);
assertTrue(sequencer.runVMAT);

function test_siochiVMATStaticModeUnchanged
% Regression guard: changing this class's ancestry must not affect static
% (non-VMAT) sequencing behavior.
p = load('photons_testData.mat');
sequencer = matRad_SequencingPhotonsSiochiLeaf(p.pln);
assertFalse(sequencer.runVMAT);

sequence = sequencer.sequence(p.resultGUI.w, p.stf);
assertTrue(isstruct(sequence.apertureInfo));
assertFalse(sequence.apertureInfo.runVMAT);
assertEqual(numel(sequence.beam), numel(p.stf));

function test_siochiVMATGatesToFMOBeamsOnly
[stf, pln] = helper_getVmatStf();
sequencer = helper_getSequencer(pln);

w = ones(sum([stf.numOfRays]), 1);
sequence = sequencer.sequence(w, stf);

fmoBeams = find(arrayfun(@(s) s.arc.isFMOBeam, stf));
assertTrue(~isempty(fmoBeams));

for i = 1:numel(stf)
    assertTrue(isfield(sequence.beam(i), 'bixelIx'));
end

function test_siochiVMATApertureInfoShape
[stf, pln] = helper_getVmatStf();
sequencer = helper_getSequencer(pln);

w = ones(sum([stf.numOfRays]), 1);
sequence = sequencer.sequence(w, stf);

apInfo = sequence.apertureInfo;
assertTrue(apInfo.runVMAT);
assertTrue(isfield(apInfo, 'arc'));
assertEqual(numel(apInfo.arc.beam), numel(stf));
assertTrue(isfield(apInfo, 'bixelWeights'));
assertEqual(numel(sequence.w), sum([stf.numOfRays]));

% every DAO beam has exactly one shape with a single set of leaf positions
for i = 1:numel(stf)
    if apInfo.arc.beam(i).isDAOBeam
        assertEqual(apInfo.beam(i).numOfShapes, 1);
        assertTrue(isfield(apInfo.beam(i).shape(1), 'leftLeafPos'));
        assertTrue(isfield(apInfo.beam(i).shape(1), 'MURate'));
    end
end

function test_siochiVMATNumOfBeamChildrenRespected
[stf, pln] = helper_getVmatStf();
sequencer = helper_getSequencer(pln);

w = ones(sum([stf.numOfRays]), 1);
sequence = sequencer.sequence(w, stf);

% every FMO beam must end up with exactly numOfChildren shapes
% distributed across its children (discardApertures caps it)
for i = 1:numel(stf)
    if stf(i).arc.isFMOBeam
        assertTrue(sequence.beam(i).numOfShapes <= stf(i).arc.numOfChildren);
    end
end

function test_siochiVMATPreconditionerRuns
[stf, pln] = helper_getVmatStf();
sequencer = helper_getSequencer(pln);
sequencer.preconditioner = true;

w = ones(sum([stf.numOfRays]), 1);
sequence = sequencer.sequence(w, stf);

assertTrue(isfield(sequence.apertureInfo, 'jacobiScale'));

function test_vmatDeliveryMetrics
[stf, pln] = helper_getVmatStf();
sequencer = helper_getSequencer(pln);

w = ones(sum([stf.numOfRays]), 1);
sequence = sequencer.sequence(w, stf);

result.apertureInfo = sequence.apertureInfo;
result = matRad_calcDeliveryMetrics(result, stf);

metrics = result.deliveryMetrics;
assertTrue(isfinite(metrics.planMU) && metrics.planMU > 0);
assertTrue(isfinite(metrics.planTime) && metrics.planTime > 0);
assertTrue(all(isfinite(metrics.leafSpeed)) && all(metrics.leafSpeed >= 0));
assertEqual(numel(metrics.MURate), result.apertureInfo.totalNumOfShapes);
assertTrue(all(metrics.fracMaxLeafSpeed <= 1) && all(metrics.fracMaxLeafSpeed >= 0));

function test_vmatDeliveryMetricsContinuousAperture
% continuous aperture must be requested via pln BEFORE stf generation so
% that the stf carries the doseAngleDAO bookkeeping; the sequencer then
% picks the flag up from pln.propSeq automatically
[stf, pln] = helper_getVmatStf(true);
sequencer = helper_getSequencer(pln);
assertTrue(sequencer.continuousAperture);

w = ones(sum([stf.numOfRays]), 1);
sequence = sequencer.sequence(w, stf);

result.apertureInfo = sequence.apertureInfo;
result = matRad_calcDeliveryMetrics(result, stf);

metrics = result.deliveryMetrics;
assertTrue(isfinite(metrics.planMU) && metrics.planMU > 0);
assertTrue(all(isfinite(metrics.leafSpeed)));
assertEqual(numel(metrics.leafSpeedGantryAngle), numel(metrics.leafSpeed));
assertTrue(all(metrics.fracForward >= 0 & metrics.fracForward <= 1));
assertElementsAlmostEqual(metrics.totalFracForward + metrics.totalFracBackward, 1);

function test_leafTouchingContinuousApertureReinvocation
% leafTouching's continuous-aperture branch (taken once shape(1) already
% carries leftLeafPosInitial/Final) must read doseAngleBorders from
% apertureInfo.arc.beam, not apertureInfo.beam, which has no such field.
% That branch is never reached by the current single call site (it creates
% those fields itself, at the end of its only invocation), so this
% re-invokes it directly to exercise the branch and guard against it
% breaking again, e.g. if a future post-DAO cleanup pass starts calling it
% a second time.
[stf, pln] = helper_getVmatStf(true);
sequencer = helper_getSequencer(pln);
w = ones(sum([stf.numOfRays]), 1);
sequence = sequencer.sequence(w, stf);
apertureInfo = sequence.apertureInfo;
assertTrue(isfield(apertureInfo.beam(1).shape(1), 'leftLeafPosInitial'));

apertureInfo = matRad_leafTouching(apertureInfo);

for i = 1:numel(apertureInfo.beam)
    shape = apertureInfo.beam(i).shape(1);
    assertTrue(all(isfinite(shape.leftLeafPosInitial)) && all(isfinite(shape.rightLeafPosInitial)));
    assertTrue(all(isfinite(shape.leftLeafPosFinal)) && all(isfinite(shape.rightLeafPosFinal)));
    assertTrue(all(shape.leftLeafPosInitial <= shape.rightLeafPosInitial));
    assertTrue(all(shape.leftLeafPosFinal <= shape.rightLeafPosFinal));
end

function test_maxLeafSpeedUsesRightLeafPositions
% maxLeafSpeed's continuous-aperture branch must compute the right-leaf
% speed from the right-leaf positions, not (as previously) from the
% left-leaf positions. Zero out every left-leaf displacement and inject a
% single, isolated right-leaf-only displacement: the correct computation
% reports a positive overall max leaf speed driven by that displacement,
% while the left/right mix-up reports exactly zero (since the left-leaf
% displacement it actually reads is zero everywhere).
[stf, pln] = helper_getVmatStf(true);
sequencer = helper_getSequencer(pln);
w = ones(sum([stf.numOfRays]), 1);
sequence = sequencer.sequence(w, stf);
apertureInfo = sequence.apertureInfo;

vec = apertureInfo.apertureVector;
nLP = apertureInfo.totalNumOfLeafPairs;
injected = false;
for i = 1:numel(apertureInfo.beam)
    if isempty(apertureInfo.arc.beam(i).leafConstMask)
        continue
    end
    n = apertureInfo.beam(i).numOfActiveLeafPairs;
    if apertureInfo.arc.beam(i).isDAOBeam
        ixLI = apertureInfo.beam(i).shape(1).vectorOffset(1) + (0:n - 1);
        ixLF = apertureInfo.beam(i).shape(1).vectorOffset(2) + (0:n - 1);
    else
        ixLI = apertureInfo.beam(apertureInfo.arc.beam(i).lastDAOBeamIx).shape(1).vectorOffset(2) + (0:n - 1);
        ixLF = apertureInfo.beam(apertureInfo.arc.beam(i).nextDAOBeamIx).shape(1).vectorOffset(1) + (0:n - 1);
    end
    ixRI = ixLI + nLP;
    ixRF = ixLF + nLP;

    vec(ixLF) = vec(ixLI); % zero left-leaf displacement everywhere
    if ~injected && apertureInfo.arc.beam(i).isDAOBeam
        vec(ixRF) = vec(ixRI) + 20; % isolated right-leaf-only displacement [mm]
        injected = true;
    else
        vec(ixRF) = vec(ixRI); % zero right-leaf displacement elsewhere
    end
end
assertTrue(injected);
apertureInfo.apertureVector = vec;

apertureInfo = matRad_calcMaxLeafSpeed(apertureInfo);
assertTrue(apertureInfo.maxLeafSpeed > 0);

function test_vmatRecalcApertureChain
% fine-angle aperture recalculation: interpolate the sequenced apertures
% onto a finer arc and recompute the bixel weights via matRad_refineApertureArc
[stf, pln] = helper_getVmatStf();
sequencer = helper_getSequencer(pln);

w = ones(sum([stf.numOfRays]), 1);
sequence = sequencer.sequence(w, stf);
apertureInfo = sequence.apertureInfo;

p = load('photons_testData.mat', 'ct', 'cst');
[stfFine, apertureInfoFine] = matRad_refineApertureArc(p.ct, p.cst, pln, apertureInfo, 7.5);

assertEqual(numel(apertureInfoFine.bixelWeights), sum([stfFine.totalNumOfBixels]));
assertTrue(all(isfinite(apertureInfoFine.bixelWeights)));
assertTrue(all(apertureInfoFine.bixelWeights >= 0));
assertTrue(any(apertureInfoFine.bixelWeights > 0));

function test_discardAperturesKeepsHighestDAPAndPreservesTotal
beam.shapes = false(1, 3, 3);
beam.shapes(:, :, 1) = [1 0 0];
beam.shapes(:, :, 2) = [1 1 1];
beam.shapes(:, :, 3) = [1 1 0];
beam.shapesWeight = [1; 10; 2];
beam.numOfShapes = 3;

newBeam = matRad_PhotonSequencerAbstract.discardApertures(beam, 2);

assertEqual(newBeam.numOfShapes, 2);
assertEqual(newBeam.shapes, double(beam.shapes(:, :, [2 3])));
oldDAP = sum(arrayfun(@(i) nnz(beam.shapes(:, :, i)) * beam.shapesWeight(i), 1:3));
newDAP = sum(arrayfun(@(i) nnz(newBeam.shapes(:, :, i)) * newBeam.shapesWeight(i), 1:2));
assertElementsAlmostEqual(newDAP, oldDAP);

function test_discardAperturesOutputShapeIsIndependentOfInput
% The output should not inherit the class of a logical shape stack nor the
% orientation of a row-vector weight list -- downstream code indexes the
% weights as a column and multiplies the shapes as doubles.
beam.shapes = false(1, 3, 3);
beam.shapes(:, :, 1) = [1 0 0];
beam.shapes(:, :, 2) = [1 1 1];
beam.shapes(:, :, 3) = [1 1 0];
beam.numOfShapes = 3;

beam.shapesWeight = [1 10 2];   % row vector
newBeam = matRad_PhotonSequencerAbstract.discardApertures(beam, 2);
assertEqual(size(newBeam.shapesWeight), [2 1]);
assertTrue(isa(newBeam.shapes, 'double'));

beam.shapesWeight = [1; 10; 2]; % column vector -> same result
newBeamCol = matRad_PhotonSequencerAbstract.discardApertures(beam, 2);
assertElementsAlmostEqual(newBeamCol.shapesWeight, newBeam.shapesWeight);

% asking for more apertures than exist clamps instead of padding
newBeamAll = matRad_PhotonSequencerAbstract.discardApertures(beam, 99);
assertEqual(newBeamAll.numOfShapes, 3);
assertEqual(size(newBeamAll.shapes, 3), 3);

function test_vmatContinuousRecalcPreservesBoundaryLeaves
[stf, pln] = helper_getVmatStf(true);
sequencer = helper_getSequencer(pln);
w = ones(sum([stf.numOfRays]), 1);
apertureInfo = sequencer.sequence(w, stf).apertureInfo;

p = load('photons_testData.mat', 'ct', 'cst');
[~, apertureInfoFine] = matRad_refineApertureArc(p.ct, p.cst, pln, apertureInfo, 7.5);

firstFine = apertureInfoFine.beam(1).shape(1);
lastFine = apertureInfoFine.beam(end).shape(1);
assertTrue(all(isfinite(apertureInfoFine.bixelWeights)));
assertTrue(all(isfinite(firstFine.leftLeafPosInitial)));
assertTrue(all(isfinite(firstFine.rightLeafPosInitial)));
assertTrue(all(isfinite(lastFine.leftLeafPosFinal)));
assertTrue(all(isfinite(lastFine.rightLeafPosFinal)));
assertElementsAlmostEqual(firstFine.leftLeafPosInitial, ...
                          apertureInfo.beam(1).shape(1).leftLeafPosInitial);
assertElementsAlmostEqual(lastFine.rightLeafPosFinal, ...
                          apertureInfo.beam(end).shape(1).rightLeafPosFinal);

% Regression: the beam-centre positions used to interpolate on the old beam
% centres alone, which do not bracket the outermost new centres -- the first
% and last refined beam came back NaN. matRad_recalcApertureInfo marks its
% output as discrete (continuousAperture = false), so every consumer reads
% exactly these fields.
for i = 1:numel(apertureInfoFine.beam)
    shape = apertureInfoFine.beam(i).shape(1);
    assertTrue(all(isfinite(shape.leftLeafPos)), ...
               sprintf('beam %d has non-finite leftLeafPos', i));
    assertTrue(all(isfinite(shape.rightLeafPos)), ...
               sprintf('beam %d has non-finite rightLeafPos', i));
    assertTrue(all(shape.leftLeafPos <= shape.rightLeafPos));
end

function test_vmatRecalcPreservesDAOAnchors
% Regression: matRad_recalcApertureInfo used to populate nextDAOBeamIx from
% stf.arc.lastDAOBeamIx (the two assignments shared a right-hand side), so
% every recalculated beam believed its next DAO anchor was its previous one
% and the interpolation between DAO control points was weighted against the
% wrong pair of apertures.
[stf, pln] = helper_getVmatStf();
sequencer = helper_getSequencer(pln);
w = ones(sum([stf.numOfRays]), 1);
apertureInfo = sequencer.sequence(w, stf).apertureInfo;

p = load('photons_testData.mat', 'ct', 'cst');
[stfFine, apertureInfoFine] = matRad_refineApertureArc(p.ct, p.cst, pln, apertureInfo, 7.5);

lastIx = [apertureInfoFine.arc.beam.lastDAOBeamIx];
nextIx = [apertureInfoFine.arc.beam.nextDAOBeamIx];

% each anchor must come from its own stf field, not from the other one
% (stf is a struct array, so the nested field needs arrayfun rather than [])
assertEqual(lastIx, arrayfun(@(b) b.arc.lastDAOBeamIx, stfFine));
assertEqual(nextIx, arrayfun(@(b) b.arc.nextDAOBeamIx, stfFine));

% and the two must genuinely differ somewhere: were next copied from last,
% as in the original bug, these would be identical for every beam
assertTrue(any(nextIx ~= lastIx));
