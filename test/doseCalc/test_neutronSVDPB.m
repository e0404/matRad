function test_suite = test_neutronSVDPB
% Tests for the neutron SVD pencil beam dose engine
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

function [ct, cst, stf, pln] = helper_neutronTestCase()
pln = struct('radiationMode', 'neutrons', 'machine', 'Generic');
pln.propDoseCalc.engine = 'SVDPB';
pln.propStf.gantryAngles = 0;
pln.propStf.couchAngles = 0;
pln.propStf.bixelWidth = 10;
pln.bioModel = 'none';
pln.multScen = 'nomScen';

builder = matRad_PhantomBuilder([30 30 30], [5 5 5], 1);
builder.addBoxOAR('Body', [150 150 150], 'HU', 0);
builder.addBoxTarget('Target', [20 20 20], 'HU', 0, 'objectives', DoseObjectives.matRad_SquaredDeviation(100, 10));
[ct, cst] = builder.getctcst();

stf = matRad_generateStf(ct, cst, pln);

function test_getNeutronEngineFromPln
pln = struct('radiationMode', 'neutrons', 'machine', 'Generic');
pln.propDoseCalc.engine = 'SVDPB';
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
assertTrue(isa(engine, 'DoseEngines.matRad_NeutronPencilBeamSVDEngine'));
assertTrue(DoseEngines.matRad_NeutronPencilBeamSVDEngine.isAvailable(pln));

function test_notAvailableWithoutKernels
% Monte Carlo only neutron machines carry no kernels
pln = struct('radiationMode', 'neutrons', 'machine', 'BNCT');
assertFalse(DoseEngines.matRad_NeutronPencilBeamSVDEngine.isAvailable(pln));

% photon machines are not accepted
pln = struct('radiationMode', 'photons', 'machine', 'Generic');
assertFalse(DoseEngines.matRad_NeutronPencilBeamSVDEngine.isAvailable(pln));

function test_forwardDoseAndKermaCorrection
[ct, cst, stf, pln] = helper_neutronTestCase();
w = ones(sum([stf(:).totalNumOfBixels]), 1);

engine = DoseEngines.matRad_NeutronPencilBeamSVDEngine(pln);
resultWater = engine.calcDoseForward(ct, cst, stf, w);

assertTrue(isequal(size(resultWater.physicalDose), ct.cubeDim));
assertTrue(all(isfinite(resultWater.physicalDose(:))));
assertTrue(all(resultWater.physicalDose(:) >= 0));
assertTrue(max(resultWater.physicalDose(:)) > 0);

% The KERMA correction cube must live on the dose grid
assertTrue(isequal(size(engine.cubeKERMAcorr{1}), engine.doseGrid.dimensions));

% Bone in the beam path (KERMA factor < 1 relative to water) lowers the dose
ctBone = ct;
ctBone.cubeHU{1}(:, :, :) = 1000;
engineBone = DoseEngines.matRad_NeutronPencilBeamSVDEngine(pln);
resultBone = engineBone.calcDoseForward(ctBone, cst, stf, w);

assertTrue(all(engineBone.cubeKERMAcorr{1}(:) < 1));
assertTrue(max(resultBone.physicalDose(:)) < max(resultWater.physicalDose(:)));

function test_doseInfluenceOnCoarserDoseGrid
[ct, cst, stf, pln] = helper_neutronTestCase();
pln.propDoseCalc.doseGrid.resolution = struct('x', 7.5, 'y', 7.5, 'z', 7.5);

dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

assertTrue(isequal(dij.doseGrid.dimensions, [20 20 20]));
assertEqual(size(dij.physicalDose{1}), [dij.doseGrid.numOfVoxels, dij.totalNumOfBixels]);
assertTrue(nnz(dij.physicalDose{1}) > 0);
