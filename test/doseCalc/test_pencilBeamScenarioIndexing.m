function test_suite = test_pencilBeamScenarioIndexing
% Regression tests for #935 and #905: the pencil beam engines have to use
% the radiological depths of the right CT scenario together with the range
% shift of the right scenario, also if only a subset of the CT scenarios is
% used
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

function [ct, cst] = helper_multiCtPhantom(numOfCtScen)
% Water box whose density differs between the CT scenarios
ct = struct();
ct.cubeDim = [16 16 16];
ct.resolution = struct('x', 4, 'y', 4, 'z', 4);
ct.numOfCtScen = numOfCtScen;
huValues = linspace(-300, 300, numOfCtScen);
ct.cubeHU = arrayfun(@(hu) hu * ones(ct.cubeDim), huValues, 'UniformOutput', false);

body = (1:prod(ct.cubeDim))';
[x, y, z] = ndgrid(7:10, 7:10, 7:10);
target = sub2ind(ct.cubeDim, x(:), y(:), z(:));

cst = {0, 'Body', 'OAR', repmat({body}, 1, numOfCtScen), struct('Priority', 2), []; ...
       1, 'Target', 'TARGET', repmat({target}, 1, numOfCtScen), struct('Priority', 1), []};

function [pln, stf] = helper_plan(ct, cst, radiationMode, engine)
pln = struct('radiationMode', radiationMode, 'machine', 'Generic', 'numOfFractions', 1);
pln.propStf.gantryAngles = 0;
pln.propStf.couchAngles = 0;
pln.propStf.bixelWidth = 8;
pln.propStf.isoCenter = matRad_getIsoCenter(cst, ct, 0);
pln.propDoseCalc.engine = engine;
pln.propDoseCalc.doseGrid.resolution = ct.resolution;
pln.propDoseCalc.enableDijSampling = false; % sampling is random
pln.multScen = 'nomScen';
stf = matRad_generateStf(ct, cst, pln);

function helper_assertScenarioDoseMatchesSingleCt(radiationMode, engine)
% Every scenario of a dose calculation on all CT scenarios has to equal the
% same scenario computed on its CT scenario alone
numOfCtScen = 3;
[ct, cst] = helper_multiCtPhantom(numOfCtScen);
[pln, stf] = helper_plan(ct, cst, radiationMode, engine);

pln.multScen = matRad_WorstCaseScenarios(ct);
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

nominalDose = cell(1, numOfCtScen);
for ctScen = 1:numOfCtScen
    refModel = matRad_WorstCaseScenarios(ct);
    refModel.ctScenProb = [ctScen 1];
    pln.multScen = refModel;
    dijRef = matRad_calcDoseInfluence(ct, cst, stf, pln);

    for s = 1:refModel.totNumScen
        scenIx = refModel.scenNum2scenIx(s);
        refDose = full(dijRef.physicalDose{scenIx});
        assertTrue(any(refDose(:) > 0));
        assertElementsAlmostEqual(full(dij.physicalDose{scenIx}), refDose, 'absolute', 1e-6 * max(refDose(:)));
    end

    % scenarios of CT scenarios not in the model are not computed
    isOtherCt = true(numOfCtScen, 1);
    isOtherCt(ctScen) = false;
    otherScen = dijRef.physicalDose(isOtherCt, :, :);
    assertTrue(all(cellfun(@isempty, otherScen(:))));

    nominalDose{ctScen} = full(dij.physicalDose{refModel.sub2scenIx(ctScen, 1, 1)});
end

% the CT scenarios differ, so must the doses
for ctScen = 2:numOfCtScen
    assertTrue(max(abs(nominalDose{ctScen}(:) - nominalDose{1}(:))) > 1e-3 * max(nominalDose{1}(:)));
end

function test_photonScenarioDoseMatchesSingleCt
helper_assertScenarioDoseMatchesSingleCt('photons', 'SVDPB');

function test_protonScenarioDoseMatchesSingleCt
helper_assertScenarioDoseMatchesSingleCt('protons', 'HongPB');

function test_sparseNominalCtScenarioDose
% Nominal scenario model on a subset of the CT scenarios
[ct, cst] = helper_multiCtPhantom(3);
[pln, stf] = helper_plan(ct, cst, 'photons', 'SVDPB');

dijAll = matRad_calcDoseInfluence(ct, cst, stf, pln);

model = matRad_NominalScenario(ct);
model.ctScenProb = [3 0.5; 2 0.5];
pln.multScen = model;
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

assertTrue(isempty(dij.physicalDose{1}));
for ctScen = [2 3]
    refDose = full(dijAll.physicalDose{ctScen});
    assertElementsAlmostEqual(full(dij.physicalDose{ctScen}), refDose, 'absolute', 1e-6 * max(refDose(:)));
end

function test_setupOnlyScenariosSingleCt
% Single CT scenario without range scenarios: scenMask is a 1 x 7 row
% vector, every active scenario has to be computed
[ct, cst] = helper_multiCtPhantom(1);
[pln, stf] = helper_plan(ct, cst, 'photons', 'SVDPB');

dijNominal = matRad_calcDoseInfluence(ct, cst, stf, pln);

model = matRad_WorstCaseScenarios(ct);
model.numOfRangeGridPoints = 1;
assertEqual(size(model.scenMask), [1 model.totNumScen]);
pln.multScen = model;
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

for scenIx = model.scenNum2scenIx()'
    assertTrue(nnz(dij.physicalDose{scenIx}) > 0);
end
refDose = full(dijNominal.physicalDose{1});
assertElementsAlmostEqual(full(dij.physicalDose{model.sub2scenIx(1, 1, 1)}), refDose, 'absolute', 1e-6 * max(refDose(:)));
