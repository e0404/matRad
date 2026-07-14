function test_suite = test_sequencerBase
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

function pln = helper_getPln(radMode)
p = load([radMode '_testData.mat'], 'pln');
pln = p.pln;

function test_getAvailableSequencers_photons
pln = helper_getPln('photons');
classList = matRad_SequencerBase.getAvailableSequencers(pln);
shortNames = {classList.shortName};
% all three photon leaf sequencers are available, the particle one is not
assertTrue(ismember('siochi', shortNames));
assertTrue(ismember('xia', shortNames));
assertTrue(ismember('engel', shortNames));
assertFalse(ismember('IMPT', shortNames));

function test_getAvailableSequencers_protons
pln = helper_getPln('protons');
classList = matRad_SequencerBase.getAvailableSequencers(pln);
shortNames = {classList.shortName};
% for particles only the particle spot sequencer is available
assertTrue(ismember('IMPT', shortNames));
assertFalse(ismember('siochi', shortNames));

function test_getSequencerFromPln_explicitSelection
pln = helper_getPln('photons');
pln.propSeq.sequencer = 'xia';
seq = matRad_SequencerBase.getSequencerFromPln(pln);
assertTrue(isa(seq, 'matRad_SequencingPhotonsXiaLeaf'));
assertEqual(seq.shortName, 'xia');

function test_getSequencerFromPln_invalidFallsBackToDefault
pln = helper_getPln('photons');
pln.propSeq.sequencer = 'thisSequencerDoesNotExist';
% warnDefault = false to keep the test output clean
seq = matRad_SequencerBase.getSequencerFromPln(pln, false);
% default photon sequencer is siochi (see MatRad_Config defaults)
assertEqual(seq.shortName, 'siochi');

function test_getSequencerFromPln_objectPassthrough
pln = helper_getPln('photons');
sequencer = matRad_SequencingPhotonsEngelLeaf(pln);
pln.propSeq = sequencer;
seq = matRad_SequencerBase.getSequencerFromPln(pln);
% a sequencer object stored in pln.propSeq must be returned as-is
assertTrue(seq == sequencer);

function test_getSequencerFromPln_particleDefault
pln = helper_getPln('protons');
seq = matRad_SequencerBase.getSequencerFromPln(pln, false);
assertTrue(isa(seq, 'matRad_ParticleSequencer'));
assertEqual(seq.shortName, 'IMPT');

function test_assignPropertiesFromPln_setsSequencingLevel
pln = helper_getPln('photons');
pln.propSeq.sequencer = 'siochi';
pln.propSeq.sequencingLevel = 7;
seq = matRad_SequencerBase.getSequencerFromPln(pln);
% properties given in pln.propSeq are transferred onto the sequencer
assertEqual(seq.sequencingLevel, 7);
