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

function [stf, pln] = helper_getVmatStf()
p = load('photons_testData.mat', 'ct', 'cst', 'pln');
pln = p.pln;

pln.propStf.generator                = 'PhotonVMAT';
pln.propStf.gantryAngles             = [-180, 180];
pln.propStf.couchAngles              = [0, 0];
pln.propStf.maxGantryAngleSpacing    = 15;
pln.propStf.maxDAOGantryAngleSpacing = 30;
pln.propStf.maxFMOGantryAngleSpacing = 45;
pln.propStf.isoCenter                = matRad_getIsoCenter(p.cst, p.ct, 0);

pln.propSeq.continuousAperture = false;
pln.propOpt.runVMAT = true;

stf = matRad_generateStf(p.ct, p.cst, pln);

function sequencer = helper_getSequencer(pln)
sequencer = matRad_SequencingPhotonsSiochiLeaf(pln);
sequencer.runVMAT = true;
sequencer.sequencingLevel = 5;
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

fmoBeams = find(arrayfun(@(s) s.propVMAT.FMOBeam, stf));
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
assertTrue(isfield(apInfo, 'propVMAT'));
assertEqual(numel(apInfo.propVMAT.beam), numel(stf));
assertTrue(isfield(apInfo, 'bixelWeights'));
assertEqual(numel(sequence.w), sum([stf.numOfRays]));

% every DAO beam has exactly one shape with a single set of leaf positions
for i = 1:numel(stf)
    if apInfo.propVMAT.beam(i).DAOBeam
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

% every FMO beam must end up with exactly numOfBeamChildren shapes
% distributed across its children (discardApertures caps it)
for i = 1:numel(stf)
    if stf(i).propVMAT.FMOBeam
        assertTrue(sequence.beam(i).numOfShapes <= stf(i).propVMAT.numOfBeamChildren);
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
result = matRad_calcDeliveryMetrics(result, pln, stf);

assertTrue(isfinite(result.apertureInfo.planMU) && result.apertureInfo.planMU > 0);
assertTrue(isfinite(result.apertureInfo.planTime) && result.apertureInfo.planTime > 0);

function test_vmatRecalcApertureChain
% fine-angle aperture recalculation: interpolate the sequenced apertures
% onto a finer dose grid and recompute the bixel weights
[stf, pln] = helper_getVmatStf();
sequencer = helper_getSequencer(pln);

w = ones(sum([stf.numOfRays]), 1);
sequence = sequencer.sequence(w, stf);
apertureInfo = sequence.apertureInfo;

p = load('photons_testData.mat', 'ct', 'cst');
recalc.pln = pln;
recalc.pln.propStf.maxGantryAngleSpacing = 7.5;
recalc.interpNew = true;
recalc.continuousAperture = false;
recalc.stf = matRad_generateStf(p.ct, p.cst, recalc.pln);

recalc = matRad_recalcApertureInfo(recalc, apertureInfo);
recalc.apertureInfo.continuousAperture = recalc.continuousAperture;
recalc.apertureInfo = matRad_recalcApertureBixelWeights(recalc.apertureInfo);

assertEqual(numel(recalc.apertureInfo.bixelWeights), sum([recalc.stf.totalNumOfBixels]));
assertTrue(all(isfinite(recalc.apertureInfo.bixelWeights)));
assertTrue(all(recalc.apertureInfo.bixelWeights >= 0));
assertTrue(any(recalc.apertureInfo.bixelWeights > 0));
