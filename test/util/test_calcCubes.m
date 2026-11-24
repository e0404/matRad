function test_suite = test_calcCubes
% Tests for matRad_calcCubes, in particular the aggregation of uncertainty
% quantities (standard deviation and variance) over bixels and beams.
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

function dij = helper_createDij(numOfBeams, beamNum, numOfVoxels)
% two bixels per beam hitting a grid of numOfVoxels voxels
dij.numOfBeams = numOfBeams;
dij.beamNum = beamNum(:);
dij.doseGrid.dimensions = [numOfVoxels 1 1];
numOfBixels = numel(beamNum);
dij.physicalDose = {ones(numOfVoxels, numOfBixels)};
dij.physicalDose_std = {0.1 * ones(numOfVoxels, numOfBixels)};
dij.physicalDose_MCvar = {0.01 * ones(numOfVoxels, numOfBixels)};

function test_calcCubesPhysicalDoseIsLinear
dij = helper_createDij(1, [1; 1], 1);
resultGUI = matRad_calcCubes([1; 2], dij);
assertElementsAlmostEqual(resultGUI.physicalDose, 3, 'absolute', 1e-12);
assertElementsAlmostEqual(resultGUI.physicalDose_beam1, 3, 'absolute', 1e-12);

function test_calcCubesStdIsAddedQuadraticallyWithWeights
dij = helper_createDij(1, [1; 1], 1);
resultGUI = matRad_calcCubes([1; 2], dij);
expected = sqrt(1^2 + 2^2) * 0.1;
assertElementsAlmostEqual(resultGUI.physicalDose_std, expected, 'absolute', 1e-12);
assertElementsAlmostEqual(resultGUI.physicalDose_std_beam1, expected, 'absolute', 1e-12);

function test_calcCubesBatchStdIsAddedQuadraticallyWithWeights
dij = helper_createDij(1, [1; 1], 1);
dij.physicalDose_batchStd = dij.physicalDose_std;
resultGUI = matRad_calcCubes([1; 2], dij);
assertElementsAlmostEqual(resultGUI.physicalDose_batchStd, sqrt(5) * 0.1, 'absolute', 1e-12);

function test_calcCubesVarianceIsWeightedWithSquaredWeights
dij = helper_createDij(1, [1; 1], 1);
resultGUI = matRad_calcCubes([1; 2], dij);
expected = (1^2 + 2^2) * 0.01;
assertElementsAlmostEqual(resultGUI.physicalDose_MCvar, expected, 'absolute', 1e-12);
assertElementsAlmostEqual(resultGUI.physicalDose_MCvar_beam1, expected, 'absolute', 1e-12);

function test_calcCubesUncertaintyPerBeam
% two beams with two bixels each
dij = helper_createDij(2, [1; 1; 2; 2], 1);
w = [1; 2; 3; 4];
resultGUI = matRad_calcCubes(w, dij);

assertElementsAlmostEqual(resultGUI.physicalDose_std_beam1, sqrt(1 + 4) * 0.1, 'absolute', 1e-12);
assertElementsAlmostEqual(resultGUI.physicalDose_std_beam2, sqrt(9 + 16) * 0.1, 'absolute', 1e-12);
assertElementsAlmostEqual(resultGUI.physicalDose_std, sqrt(1 + 4 + 9 + 16) * 0.1, 'absolute', 1e-12);

assertElementsAlmostEqual(resultGUI.physicalDose_MCvar_beam1, (1 + 4) * 0.01, 'absolute', 1e-12);
assertElementsAlmostEqual(resultGUI.physicalDose_MCvar_beam2, (9 + 16) * 0.01, 'absolute', 1e-12);
assertElementsAlmostEqual(resultGUI.physicalDose_MCvar, (1 + 4 + 9 + 16) * 0.01, 'absolute', 1e-12);

function test_calcCubesUncertaintyDoseToWaterAndSparse
dij = helper_createDij(1, [1; 1], 2);
dij.doseToWater = {sparse(dij.physicalDose{1})};
dij.doseToWater_std = {sparse(dij.physicalDose_std{1})};
dij.doseToWater_MCvar = {sparse(dij.physicalDose_MCvar{1})};
resultGUI = matRad_calcCubes([1; 2], dij);
assertEqual(size(resultGUI.doseToWater_std), [2 1]);
assertElementsAlmostEqual(resultGUI.doseToWater, [3; 3], 'absolute', 1e-12);
assertElementsAlmostEqual(resultGUI.doseToWater_std, sqrt(5) * 0.1 * [1; 1], 'absolute', 1e-12);
assertElementsAlmostEqual(resultGUI.doseToWater_MCvar, 0.05 * [1; 1], 'absolute', 1e-12);

function test_calcCubesUncertaintyNaNsAreSetToZero
dij = helper_createDij(1, [1; 1], 2);
dij.physicalDose_std{1}(1, 1) = NaN;
dij.physicalDose_MCvar{1}(2, 1) = NaN;
resultGUI = matRad_calcCubes([1; 2], dij);
assertFalse(any(isnan(resultGUI.physicalDose_std(:))));
assertFalse(any(isnan(resultGUI.physicalDose_MCvar(:))));

function test_calcCubesUncertaintyForUnitWeights
% forward calculations of Monte Carlo engines call matRad_calcCubes with ones;
% several bixel columns still add up in quadrature
dij = helper_createDij(1, [1; 1], 1);
resultGUI = matRad_calcCubes(ones(2, 1), dij);
assertElementsAlmostEqual(resultGUI.physicalDose_std, sqrt(2) * 0.1, 'absolute', 1e-12);
assertElementsAlmostEqual(resultGUI.physicalDose_MCvar, 2 * 0.01, 'absolute', 1e-12);
