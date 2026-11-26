function test_suite = test_doseGrid

test_functions = localfunctions();

initTestSuite;

function [ct, cst, pln] = helper_waterBoxSetup(doseResolution)
% Water box whose extent is not a multiple of the dose grid resolution
ct = struct();
ct.cubeDim = [22 23 21];
ct.resolution = struct('x', 3, 'y', 3, 'z', 2.5);
ct.cubeHU = {zeros(ct.cubeDim)};
ct.numOfCtScen = 1;
ct = matRad_getWorldAxes(ct);

targetIx = sub2ind(ct.cubeDim, 11, 12, 11);
cst = {0, 'Body', 'OAR', {(1:prod(ct.cubeDim))'}, struct('Priority', 2), []; ...
       1, 'Target', 'TARGET', {targetIx}, struct('Priority', 1), []};

pln = struct('radiationMode', 'photons', 'machine', 'Generic', 'numOfFractions', 1);
pln.propStf.gantryAngles = 0;
pln.propStf.couchAngles = 0;
pln.propStf.bixelWidth = 5;
pln.propStf.isoCenter = matRad_getIsoCenter(cst, ct, 0);
pln.propDoseCalc.engine = 'SVDPB';
pln.propDoseCalc.doseGrid.resolution = doseResolution;
pln.multScen = 'nomScen';

function test_doseGridIsCenteredAndCoversCt
[ct, cst, pln] = helper_waterBoxSetup(struct('x', 5, 'y', 5, 'z', 5));
stf = matRad_generateStf(ct, cst, pln);
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

grids = {dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z};
ctAxes = {ct.x, ct.y, ct.z};
resDose = [5 5 5];
resCt = [ct.resolution.x ct.resolution.y ct.resolution.z];
for i = 1:3
    assertElementsAlmostEqual((grids{i}(1) + grids{i}(end)) / 2, (ctAxes{i}(1) + ctAxes{i}(end)) / 2, 'absolute', 1e-9);
    assertTrue(grids{i}(1) - resDose(i) / 2 <= ctAxes{i}(1) - resCt(i) / 2 + 1e-9);
    assertTrue(grids{i}(end) + resDose(i) / 2 >= ctAxes{i}(end) + resCt(i) / 2 - 1e-9);
    assertElementsAlmostEqual(diff(grids{i}), resDose(i) * ones(1, numel(grids{i}) - 1), 'absolute', 1e-9);
end

assertEqual(dij.doseGrid.dimensions, [numel(dij.doseGrid.y) numel(dij.doseGrid.x) numel(dij.doseGrid.z)]);
assertEqual(dij.doseGrid.numOfVoxels, prod(dij.doseGrid.dimensions));
assertEqual(size(dij.physicalDose{1}, 1), dij.doseGrid.numOfVoxels);

function test_equalResolutionKeepsCtGrid
[ct, cst, pln] = helper_waterBoxSetup(struct('x', 3, 'y', 3, 'z', 2.5));
stf = matRad_generateStf(ct, cst, pln);
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

assertTrue(isequal(dij.doseGrid.x, ct.x));
assertTrue(isequal(dij.doseGrid.y, ct.y));
assertTrue(isequal(dij.doseGrid.z, ct.z));
assertEqual(dij.doseGrid.dimensions, ct.cubeDim);

function test_userProvidedDoseGridIsKept
[ct, cst, pln] = helper_waterBoxSetup(struct('x', 5, 'y', 5, 'z', 5));
stf = matRad_generateStf(ct, cst, pln);
pln.propDoseCalc.doseGrid.x = ct.x(1:2:end);
pln.propDoseCalc.doseGrid.y = ct.y(1:2:end);
pln.propDoseCalc.doseGrid.z = ct.z(1:2:end);
pln.propDoseCalc.doseGrid.resolution = struct('x', 6, 'y', 6, 'z', 5);
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

assertTrue(isequal(dij.doseGrid.x, ct.x(1:2:end)));
assertTrue(isequal(dij.doseGrid.y, ct.y(1:2:end)));
assertTrue(isequal(dij.doseGrid.z, ct.z(1:2:end)));

function test_resampledDoseGridForwardCalcIsFinite
[ct, cst, pln] = helper_waterBoxSetup(struct('x', 5, 'y', 5, 'z', 5));
stf = matRad_generateStf(ct, cst, pln);
resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, ones(sum([stf(:).totalNumOfBixels]), 1));

dose = resultGUI.physicalDose;
assertEqual(size(dose), ct.cubeDim);
assertTrue(all(isfinite(dose(:))));
assertTrue(max(dose(:)) > 0);

% The result on the coarse grid has to agree with the one on the CT grid
pln.propDoseCalc.doseGrid.resolution = ct.resolution;
resultCtGrid = matRad_calcDoseForward(ct, cst, stf, pln, ones(sum([stf(:).totalNumOfBixels]), 1));
assertElementsAlmostEqual(max(dose(:)), max(resultCtGrid.physicalDose(:)), 'relative', 0.1);

% Border voxels of the CT must not be zeroed or garbage by the resampling
% (the dose grid reaches beyond the CT centers on the high end)
border = dose(:, [1 end], :);
borderRef = resultCtGrid.physicalDose(:, [1 end], :);
assertTrue(all(isfinite(border(:))));
assertElementsAlmostEqual(border, borderRef, 'absolute', 0.15 * max(resultCtGrid.physicalDose(:)));

function test_ctOnDoseGridHasNoNaNs
[ct, cst, pln] = helper_waterBoxSetup(struct('x', 5, 'y', 5, 'z', 5));
stf = matRad_generateStf(ct, cst, pln);
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

ct = matRad_calcWaterEqD(ct, stf(1).radiationMode);
ctR = matRad_resampleCTtoGrid(ct, dij);
assertTrue(all(isfinite(ctR.cubeHU{1}(:))));
assertTrue(all(isfinite(ctR.cube{1}(:))));
