function test_suite = test_particleSequencing
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

function [resultGUI, stf, dij, pln] = helper_getTestData()
p = load('protons_testData.mat');
pln = p.pln;
resultGUI = p.resultGUI;
stf = p.stf;
dij = p.dij;

function test_particleSequencingKeepsDoseWithDij
% Regression test: for particle plans the sequencer only derives the spot
% delivery order/timing and does not modify the fluence. Passing dij must not
% trigger a dose recomputation (which used to crash on the per-beam struct
% array) and must leave the original resultGUI dose untouched.
[resultGUI, stf, dij, pln] = helper_getTestData();
fn_old = fieldnames(resultGUI);

resultGUI_sequenced = matRad_sequencing(resultGUI, stf, pln, dij);

% All original fields are preserved and unchanged (fluence/dose not modified)
fn_new = fieldnames(resultGUI_sequenced);
for i = 1:numel(fn_old)
    assertTrue(any(strcmp(fn_old{i}, fn_new)));
    assertElementsAlmostEqual(resultGUI.(fn_old{i}), resultGUI_sequenced.(fn_old{i}));
end

% Particle sequencing returns a per-beam struct array without aperture info
seq = resultGUI_sequenced.sequencing;
assertTrue(isstruct(seq));
assertEqual(numel(seq), numel(stf));
assertFalse(isfield(seq, 'apertureInfo'));
for i = 1:numel(seq)
    assertTrue(isvector(seq(i).orderToSTF));
    assertTrue(isvector(seq(i).orderToSS));
    assertTrue(isvector(seq(i).time));
    assertTrue(isvector(seq(i).w));
end

function test_particleSequencingWithoutDij
% Sequencing without dij must also work and only attach the sequencing info.
[resultGUI, stf, dij, pln] = helper_getTestData();

resultGUI_sequenced = matRad_sequencing(resultGUI, stf, pln);

assertTrue(isstruct(resultGUI_sequenced.sequencing));
assertEqual(numel(resultGUI_sequenced.sequencing), numel(stf));

function test_particleMakePhaseMatrix
% Covers the spot-time based phase matrix generation used in 4D dose calc.
[resultGUI, stf, dij, pln] = helper_getTestData();

numOfPhases  = 4;
motionPeriod = 5; % [s]

sequencer = matRad_ParticleScanningSequencerSpill();
sequence  = sequencer.sequence(resultGUI.w, stf);
sequence  = sequencer.makePhaseMatrix(sequence, numOfPhases, motionPeriod);

for i = 1:numel(sequence)
    % phase matrix has one row per spot and one column per phase
    assertEqual(size(sequence(i).phaseMatrix, 1), numel(sequence(i).time));
    assertEqual(size(sequence(i).phaseMatrix, 2), numOfPhases);
    % every spot is assigned to exactly one phase
    assertEqual(numel(sequence(i).phaseNum), numel(sequence(i).time));
end
