function test_suite = test_MCNPEngine
% Tests for the MCNP neutron dose engine (run file export only, MCNP itself
% is not required)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_functions = localfunctions();

initTestSuite;

function [ct, cst, stf, pln] = helper_mcnpTestCase()
pln = struct('radiationMode', 'neutrons', 'machine', 'BNCT');
pln.propDoseCalc.engine = 'MCNP';
pln.propDoseCalc.doseGrid.resolution = struct('x', 5, 'y', 5, 'z', 5);
pln.propStf.gantryAngles = 0;
pln.propStf.couchAngles = 0;
pln.propStf.bixelWidth = 20;
pln.bioModel = 'none';
pln.multScen = 'nomScen';

builder = matRad_PhantomBuilder([20 20 20], [5 5 5], 1);
builder.addBoxOAR('Body', [100 100 100], 'HU', 0);
builder.addBoxTarget('Target', [20 20 20], 'HU', 0, 'objectives', DoseObjectives.matRad_SquaredDeviation(100, 10));
[ct, cst] = builder.getctcst();

stf = matRad_generateStf(ct, cst, pln);

function test_getMCNPEngineFromPln
pln = struct('radiationMode', 'neutrons', 'machine', 'BNCT');
pln.propDoseCalc.engine = 'MCNP';
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
assertTrue(isa(engine, 'DoseEngines.matRad_NeutronMCNPEngine'));
assertTrue(DoseEngines.matRad_NeutronMCNPEngine.isAvailable(pln));
assertTrue(strcmp(engine.externalCalculation, 'write'));

% constant RBE is taken from the biological model
pln.bioModel = matRad_ConstantRBE();
pln.bioModel.RBE = 4;
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
assertEqual(engine.constantRBE, 4);

function test_writeRunFiles
% The tissue segmentation needs image processing and statistics functions
% that may be unavailable (e.g. Octave without the respective packages)
requiredFuncs = {'imboxfilt3', 'knnsearch', 'bwconncomp'};
if any(cellfun(@(f) isempty(which(f)), requiredFuncs))
    moxunit_throw_test_skipped_exception('Image processing/statistics functions required for the tissue segmentation not available');
end

[ct, cst, stf, pln] = helper_mcnpTestCase();

engine = DoseEngines.matRad_NeutronMCNPEngine(pln);
engine.externalCalculation = 'write';
engine.useDICOMinfoRescale = false;
engine.workingDir = helper_temporaryFolder('testMCNP', true);

dij = engine.calcDoseInfluence(ct, cst, stf);

runDir = fullfile(engine.workingDir, 'runfiles');
assertTrue(isfolder(runDir));
runFiles = dir(fullfile(runDir, 'MCNPrunfile_*bixel'));
assertEqual(numel(runFiles), dij.totalNumOfBixels);
assertTrue(numel(dir(fullfile(runDir, 'runAll.*'))) == 1);

% empty dose container is returned in write mode
assertEqual(size(dij.physicalDose{1}), [dij.doseGrid.numOfVoxels, dij.totalNumOfBixels]);
assertEqual(nnz(dij.physicalDose{1}), 0);

function test_missingBodyStructure
[ct, cst, stf, pln] = helper_mcnpTestCase();
cst{1, 2} = 'Patient';

engine = DoseEngines.matRad_NeutronMCNPEngine(pln);
engine.externalCalculation = 'write';
engine.useDICOMinfoRescale = false;
engine.workingDir = helper_temporaryFolder('testMCNP', true);

assertExceptionThrown(@() engine.calcDoseInfluence(ct, cst, stf));
