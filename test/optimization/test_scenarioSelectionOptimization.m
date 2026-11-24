function test_suite = test_scenarioSelectionOptimization
% Regression tests for the selection of scenarios and their probabilities
% in matRad_fluenceOptimization with multiple and sparse CT scenarios
% (#935, #905). The dose influence is replaced by a single voxel and a
% single bixel with a different dose per scenario, so that the optimal
% weight is known analytically.
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

function [dij, cst, pln, scenDose] = helper_singleVoxelProblem(model, robustness)
% Dose influence on all scenarios of the scenario model with a single
% influencing bixel on a single target voxel
ct = struct();
ct.cubeDim = [10 10 10];
ct.resolution = struct('x', 5, 'y', 5, 'z', 5);
ct.numOfCtScen = model.numOfAvailableCtScen;
ct.cubeHU = repmat({zeros(ct.cubeDim)}, 1, ct.numOfCtScen);

targetVoxel = sub2ind(ct.cubeDim, 5, 5, 5);
[x, y, z] = ndgrid(4:6, 4:6, 4:6);
targetDoseCalc = sub2ind(ct.cubeDim, x(:), y(:), z(:));
cst = {0, 'Target', 'TARGET', repmat({targetDoseCalc}, 1, ct.numOfCtScen), struct('Priority', 1, 'alphaX', 0.1, 'betaX', 0.05), []};

pln = struct('radiationMode', 'photons', 'machine', 'Generic', 'numOfFractions', 1);
pln.propStf.gantryAngles = 0;
pln.propStf.couchAngles = 0;
pln.propStf.bixelWidth = 10;
pln.propStf.isoCenter = matRad_getIsoCenter(cst, ct, 0);
pln.propDoseCalc.engine = 'SVDPB';
pln.propDoseCalc.doseGrid.resolution = ct.resolution;
pln.propDoseCalc.enableDijSampling = false; % sampling is random
pln.multScen = model;
pln.bioModel = 'none';
pln.propOpt.quantityOpt = 'physicalDose';
pln.propOpt.maxIter = 200;

stf = matRad_generateStf(ct, cst, pln);
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

% Replace the dose influence: dose of the scenario depends on its full
% scenario index
numOfVoxels = dij.doseGrid.numOfVoxels;
scenDose = zeros(size(dij.physicalDose));
for scenIx = model.scenNum2scenIx()'
    scenDose(scenIx) = 1 + 0.25 * scenIx;
    dij.physicalDose{scenIx} = sparse(targetVoxel, 1, scenDose(scenIx), numOfVoxels, dij.totalNumOfBixels);
end

objective = DoseObjectives.matRad_SquaredDeviation(100, 2);
objective.robustness = robustness;
cst{1, 4} = repmat({targetVoxel}, 1, ct.numOfCtScen);
cst{1, 6} = {objective};

function test_stochasticOptimizationUsesMatchingProbabilities
% STOCH: f = sum_s p_s * (D_s * w - dRef)^2, minimized by
% w = dRef * sum_s p_s * D_s / sum_s p_s * D_s^2
ct.numOfCtScen = 3;
model = matRad_WorstCaseScenarios(ct);
model.ctScenProb = [1 0.2; 2 0.3; 3 0.5];
[dij, cst, pln, scenDose] = helper_singleVoxelProblem(model, 'STOCH');
pln.propOpt.scen4D = 'all';

resultGUI = matRad_fluenceOptimization(dij, cst, pln);

p = model.scenProb;
D = scenDose(model.scenNum2scenIx());
wExpected = 2 * sum(p .* D) / sum(p .* D.^2);
assertElementsAlmostEqual(resultGUI.w(1), wExpected, 'relative', 1e-4);

function test_stochasticOptimizationDefaultsToFirstCtScenario
% Without pln.propOpt.scen4D, only the scenarios of the first CT scenario
% are used
ct.numOfCtScen = 3;
model = matRad_WorstCaseScenarios(ct);
[dij, cst, pln, scenDose] = helper_singleVoxelProblem(model, 'STOCH');

resultGUI = matRad_fluenceOptimization(dij, cst, pln);

isFirstCt = model.linearMask(:, 1) == 1;
p = model.scenProb(isFirstCt);
D = scenDose(model.scenNum2scenIx(find(isFirstCt)));
wExpected = 2 * sum(p .* D) / sum(p .* D.^2);
assertElementsAlmostEqual(resultGUI.w(1), wExpected, 'relative', 1e-4);

function test_nominalOptimizationOnSparseCtScenarios
% Only CT scenario 2 is used, which is then also the default for
% optimization: w = dRef / D_2
ct.numOfCtScen = 3;
model = matRad_NominalScenario(ct);
model.ctScenProb = [2 1];
[dij, cst, pln, scenDose] = helper_singleVoxelProblem(model, 'none');

resultGUI = matRad_fluenceOptimization(dij, cst, pln);

assertElementsAlmostEqual(resultGUI.w(1), 2 / scenDose(2), 'relative', 1e-4);
assertElementsAlmostEqual(max(resultGUI.physicalDose(:)), 2, 'relative', 1e-4);

function test_nominalOptimizationOnSelectedCtScenarios
% Nominal objective summed over the selected CT scenarios 2 and 3:
% w = dRef * (D_2 + D_3) / (D_2^2 + D_3^2)
ct.numOfCtScen = 3;
model = matRad_WorstCaseScenarios(ct);
[dij, cst, pln, scenDose] = helper_singleVoxelProblem(model, 'none');
pln.propOpt.scen4D = [2 3];

resultGUI = matRad_fluenceOptimization(dij, cst, pln);

D = scenDose([2 3]);
wExpected = 2 * sum(D) / sum(D.^2);
assertElementsAlmostEqual(resultGUI.w(1), wExpected, 'relative', 1e-4);
