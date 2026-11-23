function test_suite = test_electronDensitiesToHU
% Tests for matRad_electronDensitiesToHU, in particular for the HLUT lookup
% that is required if the ct does not provide an HLUT (issue #941)
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

function ct = helper_createWaterCt()
% water phantom with relative electron density 1 everywhere and no HLUT
builder = matRad_PhantomBuilder([20 20 20], [5 5 5], 1);
builder.addBoxOAR('Water', [16 16 16], 'HU', 0);
[ct, ~] = builder.getctcst();
ct = rmfield(ct, 'cubeHU');
ct.cube = {ones(ct.cubeDim)};

function hlut = helper_bijectiveHlut()
% the matRad default HLUT has a duplicate density at its lower end, which
% makes it non-bijective. Remove the first row to obtain a bijective one.
hlut = matRad_readHLUT([MatRad_Config.instance().matRadSrcRoot filesep 'hluts' filesep 'matRad_default.hlut']);
hlut = hlut(2:end, :);

function test_noRadiationModeDoesNotFail
% regression test for #941: matRad_loadHLUT was called without radiationMode
ct = helper_createWaterCt();
ctOut = matRad_electronDensitiesToHU(ct);
assertTrue(isstruct(ctOut));
assertEqual(ctOut.cube, ct.cube);
% the default photon HLUT maps water (density 1) to 0 HU
assertTrue(isfield(ctOut, 'cubeHU'));
assertElementsAlmostEqual(ctOut.cubeHU{1}, zeros(ct.cubeDim), 'absolute', 1e-10);

function test_defaultHlutAirPlateauUsesLowestHU
% the default HLUT maps -1024 and -999 HU to the density 0.001
ct = helper_createWaterCt();
ct.cube{1}(1) = 0.001;
ctOut = matRad_electronDensitiesToHU(ct, 'photons');
assertElementsAlmostEqual(ctOut.cubeHU{1}(1), -1024, 'absolute', 1e-10);
assertElementsAlmostEqual(ctOut.cubeHU{1}(2), 0, 'absolute', 1e-10);

function test_noRadiationModeEqualsPhotons
ct = helper_createWaterCt();
assertEqual(matRad_electronDensitiesToHU(ct), matRad_electronDensitiesToHU(ct, 'photons'));
assertEqual(matRad_electronDensitiesToHU(ct), matRad_electronDensitiesToHU(ct, []));

function test_explicitRadiationModes
ct = helper_createWaterCt();
modes = {'protons', 'carbon', struct('radiationMode', 'protons')};
for i = 1:numel(modes)
    ctOut = matRad_electronDensitiesToHU(ct, modes{i});
    assertTrue(isstruct(ctOut));
    assertEqual(ctOut.cube, ct.cube);
end

function test_invalidRadiationModeThrows
ct = helper_createWaterCt();
assertExceptionThrown(@() matRad_electronDensitiesToHU(ct, 42));

function test_waterGivesZeroHUWithCtHlut
% an hlut in the ct is used, radiationMode is not needed
ct = helper_createWaterCt();
ct.hlut = helper_bijectiveHlut();
ctOut = matRad_electronDensitiesToHU(ct);
assertTrue(isfield(ctOut, 'cubeHU'));
assertEqual(size(ctOut.cubeHU{1}), reshape(ct.cubeDim, 1, []));
assertElementsAlmostEqual(ctOut.cubeHU{1}, zeros(ct.cubeDim), 'absolute', 1e-10);

function test_roundTripHUToDensityAndBack
% matRad_calcWaterEqD stores the used HLUT in the ct; make it bijective and
% convert back
ct = helper_createWaterCt();
ct = rmfield(ct, 'cube');
huRange = linspace(-990, 2900, 20 * 20 * 20);
ct.cubeHU = {reshape(huRange, ct.cubeDim)};
ct = matRad_calcWaterEqD(ct, 'photons');
assertTrue(isfield(ct, 'hlut'));
huOriginal = ct.cubeHU{1};

ct.hlut = helper_bijectiveHlut();
ct = rmfield(ct, 'cubeHU');
ct = matRad_electronDensitiesToHU(ct);
assertTrue(isfield(ct, 'cubeHU'));
assertElementsAlmostEqual(ct.cubeHU{1}, huOriginal, 'absolute', 1e-6);

function test_roundTripMultipleScenarios
ct = helper_createWaterCt();
ct.numOfCtScen = 2;
ct.cube = {0.5 * ones(ct.cubeDim), 1.2 * ones(ct.cubeDim)};
ct.hlut = helper_bijectiveHlut();
ctOut = matRad_electronDensitiesToHU(ct);
assertEqual(numel(ctOut.cubeHU), 2);
ctBack = matRad_calcWaterEqD(rmfield(ctOut, {'hlut', 'cube'}), 'photons');
assertElementsAlmostEqual(ctBack.cube{1}, ct.cube{1}, 'absolute', 1e-6);
assertElementsAlmostEqual(ctBack.cube{2}, ct.cube{2}, 'absolute', 1e-6);

function test_outOfRangeDensitiesAreProjected
ct = helper_createWaterCt();
ct.hlut = helper_bijectiveHlut();
ct.cube{1}(1) = 100;
ct.cube{1}(2) = -5;
ctOut = matRad_electronDensitiesToHU(ct);
assertElementsAlmostEqual(ctOut.cubeHU{1}(1), max(ct.hlut(:, 1)), 'absolute', 1e-10);
assertElementsAlmostEqual(ctOut.cubeHU{1}(2), min(ct.hlut(:, 1)), 'absolute', 1e-10);

function test_densityPlateauUsesLowestHU
% several HU values with the same density: the lowest HU value is used and
% the remaining HLUT is interpolated as usual
ct = helper_createWaterCt();
ct.hlut = [-1000 0.001; -500 0.5; 0 1; 50 1; 100 1; 1000 2];
ct.cube{1}(1) = 0.75;
ct.cube{1}(2) = 1.5;
ctOut = matRad_electronDensitiesToHU(ct);
assertElementsAlmostEqual(ctOut.cubeHU{1}(1), -250, 'absolute', 1e-10);
assertElementsAlmostEqual(ctOut.cubeHU{1}(2), 500, 'absolute', 1e-10);
assertElementsAlmostEqual(ctOut.cubeHU{1}(3), 0, 'absolute', 1e-10);

function test_unsortedHlutIsSorted
ct = helper_createWaterCt();
ct.hlut = flipud(helper_bijectiveHlut());
ctOut = matRad_electronDensitiesToHU(ct);
assertElementsAlmostEqual(ctOut.cubeHU{1}, zeros(ct.cubeDim), 'absolute', 1e-10);

function test_decreasingDensitiesSkipReconversion
ct = helper_createWaterCt();
ct.hlut = [-1000 0.001; 0 1; 100 0.9; 1000 2];
ctOut = matRad_electronDensitiesToHU(ct);
assertFalse(isfield(ctOut, 'cubeHU'));
assertEqual(ctOut.cube, ct.cube);
