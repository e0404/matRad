function test_suite = test_calculateProbabilisticQuantities
% Tests for matRad_calculateProbabilisticQuantities, in particular that the
% scenario probabilities are paired with the dose influence of the same
% scenario for multiple and sparse CT scenarios (#935)
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

function [dij, cst, pln] = helper_syntheticProblem(model)
numOfVoxels = 6;
numOfBixels = 2;
voi = [2; 3; 5];

dij = struct();
dij.doseGrid.numOfVoxels = numOfVoxels;
dij.totalNumOfBixels = numOfBixels;
dij.physicalDose = cell(size(model.scenMask));
for scenIx = model.scenNum2scenIx()'
    % dose influence depends on the full scenario index
    dij.physicalDose{scenIx} = sparse(scenIx * reshape(1:numOfVoxels * numOfBixels, numOfVoxels, numOfBixels));
end

objective = DoseObjectives.matRad_SquaredDeviation(1, 1);
objective.robustness = 'PROB';
cst = {0, 'Target', 'TARGET', repmat({voi}, 1, model.numOfAvailableCtScen), struct(), {objective}};

pln = struct();
pln.multScen = model;
pln.bioModel = struct('bioOpt', false);

function [expDose, omega] = helper_reference(dij, cst, model, ctScenIds)
% Reference with explicit pairing of scenario number and dose influence
voi = cst{1, 4}{1};
expDose = sparse(dij.doseGrid.numOfVoxels, dij.totalNumOfBixels);
omega = zeros(dij.totalNumOfBixels);
for s = 1:model.totNumScen
    if ~ismember(model.linearMask(s, 1), ctScenIds)
        continue
    end
    D = dij.physicalDose{model.sub2scenIx(model.linearMask(s, 1), model.linearMask(s, 2), model.linearMask(s, 3))};
    p = model.scenProb(s);
    expDose = expDose + p * D;
    omega = omega + (p * D(voi, :))' * (p * D(voi, :));
end
omega = omega - expDose(voi, :)' * expDose(voi, :);

function test_phaseModeMultipleCtScenarios
ct.numOfCtScen = 3;
model = matRad_WorstCaseScenarios(ct);
model.ctScenProb = [1 0.2; 2 0.3; 3 0.5];
[dij, cst, pln] = helper_syntheticProblem(model);

dij = matRad_calculateProbabilisticQuantities(dij, cst, pln, 'phase');

for ctScen = 1:3
    [expDose, omega] = helper_reference(dij, cst, model, ctScen);
    assertElementsAlmostEqual(full(dij.physicalDoseExp{ctScen}), full(expDose));
    assertElementsAlmostEqual(dij.physicalDoseOmega{1, ctScen}, omega);
end

function test_phaseModeSparseCtScenarios
ct.numOfCtScen = 3;
model = matRad_WorstCaseScenarios(ct);
model.ctScenProb = [3 0.5; 2 0.5];
[dij, cst, pln] = helper_syntheticProblem(model);

dij = matRad_calculateProbabilisticQuantities(dij, cst, pln);

assertTrue(isempty(dij.physicalDoseExp{1}));
for ctScen = [2 3]
    [expDose, omega] = helper_reference(dij, cst, model, ctScen);
    assertElementsAlmostEqual(full(dij.physicalDoseExp{ctScen}), full(expDose));
    assertElementsAlmostEqual(dij.physicalDoseOmega{1, ctScen}, omega);
end

function test_allModeMultipleCtScenarios
ct.numOfCtScen = 2;
model = matRad_RandomScenarios(ct);
model.nSamples = 4;
[dij, cst, pln] = helper_syntheticProblem(model);

dij = matRad_calculateProbabilisticQuantities(dij, cst, pln, 'all');

[expDose, omega] = helper_reference(dij, cst, model, [1 2]);
assertEqual(numel(dij.physicalDoseExp), 1);
assertElementsAlmostEqual(full(dij.physicalDoseExp{1}), full(expDose));
assertElementsAlmostEqual(dij.physicalDoseOmega{1, 1}, omega);
