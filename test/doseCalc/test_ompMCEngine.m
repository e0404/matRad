function test_suite = test_ompMCEngine
% Tests for the ompMC photon Monte Carlo dose engine. Tests are skipped if the
% compiled ompMC mex interface is not available on the current platform.
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

function [ct, cst, stf, pln] = helper_setupTwoBeamPlan(gantryAngles)
% Creates a small two-beam plan (several bixels per beam) on a coarse dose grid
if ~matRad_checkMexFileExists('omc_matrad')
    moxunit_throw_test_skipped_exception('ompMC mex interface not available');
end

load('BOXPHANTOM.mat', 'ct', 'cst');

pln = struct('radiationMode', 'photons', 'machine', 'Generic', 'numOfFractions', 1);
pln.propStf.gantryAngles = gantryAngles;
pln.propStf.couchAngles = zeros(size(gantryAngles));
pln.propStf.bixelWidth = 20;
pln.propStf.isoCenter = matRad_getIsoCenter(cst, ct, 0);
pln.propDoseCalc.engine = 'ompMC';
pln.propDoseCalc.doseGrid.resolution = struct('x', 10, 'y', 10, 'z', 10);
pln.propDoseCalc.numHistoriesPerBeamlet = 500;
pln.propDoseCalc.numHistoriesDirect = 500;

stf = matRad_generateStf(ct, cst, pln);

% Make sure the regression scenario (more than one bixel per beam) is covered
assertTrue(all([stf.totalNumOfBixels] > 1));

% The forward calculation should use the same number of histories per bixel
pln.propDoseCalc.numHistoriesDirect = sum([stf.totalNumOfBixels]) * pln.propDoseCalc.numHistoriesPerBeamlet;

function test_ompMCEngineAvailable
if ~matRad_checkMexFileExists('omc_matrad')
    moxunit_throw_test_skipped_exception('ompMC mex interface not available');
end
pln = struct('radiationMode', 'photons', 'machine', 'Generic');
assertTrue(DoseEngines.matRad_PhotonOmpMCEngine.isAvailable(pln));

function test_ompMCForwardMultipleBixelsPerBeam
[ct, cst, stf, pln] = helper_setupTwoBeamPlan([0 90]);

w = ones(sum([stf.totalNumOfBixels]), 1);
resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, w);

assertTrue(isfield(resultGUI, 'physicalDose'));
assertEqual(size(resultGUI.physicalDose), ct.cubeDim);
assertTrue(all(isfinite(resultGUI.physicalDose(:))));
assertTrue(all(resultGUI.physicalDose(:) >= 0));
assertTrue(max(resultGUI.physicalDose(:)) > 0);

% the total dose has to be the sum of the individual beam doses
doseBeamSum = zeros(size(resultGUI.physicalDose));
for beamIx = 1:numel(stf)
    fieldName = sprintf('physicalDose_beam%d', beamIx);
    assertTrue(isfield(resultGUI, fieldName));
    assertEqual(size(resultGUI.(fieldName)), size(resultGUI.physicalDose));
    assertTrue(max(resultGUI.(fieldName)(:)) > 0);
    doseBeamSum = doseBeamSum + resultGUI.(fieldName);
end
assertElementsAlmostEqual(doseBeamSum, resultGUI.physicalDose, 'relative', 1e-6, 1e-6);

function test_ompMCForwardMatchesDijBasedDose
[ct, cst, stf, pln] = helper_setupTwoBeamPlan([0 90]);

% non-uniform weights to also check the weighting of the individual bixels
nBixels = sum([stf.totalNumOfBixels]);
w = (1:nBixels)' ./ nBixels;

resultForward = matRad_calcDoseForward(ct, cst, stf, pln, w);

dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
resultDij = matRad_calcCubes(w, dij);

assertEqual(size(resultForward.physicalDose), size(resultDij.physicalDose));
for beamIx = 1:numel(stf)
    fieldName = sprintf('physicalDose_beam%d', beamIx);
    assertEqual(size(resultForward.(fieldName)), size(resultDij.(fieldName)));
end

% Monte Carlo noise is not identical between the two modes, so use a loose tolerance
totalForward = sum(resultForward.physicalDose(:));
totalDij = sum(resultDij.physicalDose(:));
assertTrue(abs(totalForward - totalDij) < 0.1 * totalDij);

function test_ompMCForwardMCvariance
[ct, cst, stf, pln] = helper_setupTwoBeamPlan([0 90]);
pln.propDoseCalc.outputMCvariance = true;

nBixels = sum([stf.totalNumOfBixels]);
w = (1:nBixels)' ./ nBixels;
resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, w);

assertTrue(isfield(resultGUI, 'physicalDose_MCvar'));
assertEqual(size(resultGUI.physicalDose_MCvar), size(resultGUI.physicalDose));
assertTrue(all(isfinite(resultGUI.physicalDose_MCvar(:))));
assertTrue(all(resultGUI.physicalDose_MCvar(:) >= 0));
assertTrue(max(resultGUI.physicalDose_MCvar(:)) > 0);

% independent beams: variance of the total dose is the sum of the beam variances
varBeamSum = zeros(size(resultGUI.physicalDose_MCvar));
for beamIx = 1:numel(stf)
    fieldName = sprintf('physicalDose_MCvar_beam%d', beamIx);
    assertTrue(isfield(resultGUI, fieldName));
    assertEqual(size(resultGUI.(fieldName)), size(resultGUI.physicalDose));
    varBeamSum = varBeamSum + resultGUI.(fieldName);
end
assertElementsAlmostEqual(varBeamSum, resultGUI.physicalDose_MCvar, 'relative', 1e-6, 1e-6);

function test_ompMCForwardSingleBeam
[ct, cst, stf, pln] = helper_setupTwoBeamPlan(0);

w = ones(sum([stf.totalNumOfBixels]), 1);
resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, w);

assertEqual(size(resultGUI.physicalDose), ct.cubeDim);
assertElementsAlmostEqual(resultGUI.physicalDose_beam1, resultGUI.physicalDose, 'relative', 1e-6, 1e-6);
