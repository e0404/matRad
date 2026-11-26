function test_suite = test_createDoseGridAxes

test_functions = localfunctions();

initTestSuite;

function ctGrid = helper_ctGrid(cubeDim, resolution)
ct = helper_createTestCt(cubeDim, resolution, 'createCoordinateArrays', true);
ctGrid = struct('x', ct.x, 'y', ct.y, 'z', ct.z, 'resolution', ct.resolution);

function test_doseGridIsCenteredOnCt
ctGrid = helper_ctGrid([22 23 21], [3 3 2.5]);
doseRes = struct('x', 5, 'y', 5, 'z', 5);
[x, y, z] = matRad_createDoseGridAxes(ctGrid, doseRes);

assertElementsAlmostEqual((x(1) + x(end)) / 2, (ctGrid.x(1) + ctGrid.x(end)) / 2, 'absolute', 1e-9);
assertElementsAlmostEqual((y(1) + y(end)) / 2, (ctGrid.y(1) + ctGrid.y(end)) / 2, 'absolute', 1e-9);
assertElementsAlmostEqual((z(1) + z(end)) / 2, (ctGrid.z(1) + ctGrid.z(end)) / 2, 'absolute', 1e-9);

function test_doseGridCoversCtAndKeepsResolution
ctGrid = helper_ctGrid([22 23 21], [3 3 2.5]);
doseRes = struct('x', 5, 'y', 4, 'z', 7);
[x, y, z] = matRad_createDoseGridAxes(ctGrid, doseRes);

axesDose = {x, y, z};
axesCt = {ctGrid.x, ctGrid.y, ctGrid.z};
resDose = [doseRes.x doseRes.y doseRes.z];
resCt = [ctGrid.resolution.x ctGrid.resolution.y ctGrid.resolution.z];
tol = 1e-9;
for i = 1:3
    assertTrue(isrow(axesDose{i}));
    assertElementsAlmostEqual(diff(axesDose{i}), resDose(i) * ones(1, numel(axesDose{i}) - 1), 'absolute', tol);
    assertTrue(axesDose{i}(1) - resDose(i) / 2 <= axesCt{i}(1) - resCt(i) / 2 + tol);
    assertTrue(axesDose{i}(end) + resDose(i) / 2 >= axesCt{i}(end) + resCt(i) / 2 - tol);
    % tight: never more voxels than necessary to cover the CT
    nCt = numel(axesCt{i});
    assertEqual(numel(axesDose{i}), ceil(nCt * resCt(i) / resDose(i) - 1e-6));
end

function test_equalResolutionReproducesCtGridExactly
ct = helper_createTestCt([22 23 21], [0.1 0.3 0.7], 'createCoordinateArrays', true);
ctGrid = struct('x', ct.x + 1 / 3, 'y', ct.y - 1 / 7, 'z', ct.z + pi, 'resolution', ct.resolution);
[x, y, z] = matRad_createDoseGridAxes(ctGrid, ct.resolution);

assertTrue(isequal(x, ctGrid.x));
assertTrue(isequal(y, ctGrid.y));
assertTrue(isequal(z, ctGrid.z));

function test_integerMultipleIsNotRoundedUp
% Ratios that are integers mathematically must not be rounded up by floating point noise
ctGrid = helper_ctGrid([100 100 100], [0.1 0.1 0.1]);
doseRes = struct('x', 0.3, 'y', 0.3, 'z', 0.3);
[x, y, z] = matRad_createDoseGridAxes(ctGrid, doseRes);
assertEqual(numel(x), 34); % ceil(100 * 0.1 / 0.3) = 34 (33.33...)

ctGrid = helper_ctGrid([30 30 30], [0.1 0.1 0.1]);
[x, y, z] = matRad_createDoseGridAxes(ctGrid, doseRes);
assertEqual(numel(x), 10); % 30 * 0.1 / 0.3 = 10.000000000000002 in floating point
assertEqual(numel(y), 10);
assertEqual(numel(z), 10);

function test_tg119LikeGrid
% TG119: 167 x 167 x 129 voxels of 3 x 3 x 2.5 mm, 5 mm dose grid
ctGrid = helper_ctGrid([167 167 129], [3 3 2.5]);
doseRes = struct('x', 5, 'y', 5, 'z', 5);
[x, y, z] = matRad_createDoseGridAxes(ctGrid, doseRes);

assertEqual(numel(x), 101);
assertEqual(numel(y), 101);
assertEqual(numel(z), ceil(129 * 2.5 / 5));
assertElementsAlmostEqual((x(1) + x(end)) / 2, (ctGrid.x(1) + ctGrid.x(end)) / 2, 'absolute', 1e-9);
assertElementsAlmostEqual((z(1) + z(end)) / 2, (ctGrid.z(1) + ctGrid.z(end)) / 2, 'absolute', 1e-9);

function test_singleVoxelAxis
ctGrid = helper_ctGrid([1 5 5], [3 3 3]);
doseRes = struct('x', 5, 'y', 5, 'z', 5);
[x, y, z] = matRad_createDoseGridAxes(ctGrid, doseRes);
assertEqual(numel(y), 1); % cubeDim(1) = 1 -> a single voxel along y
assertElementsAlmostEqual(y, ctGrid.y(1), 'absolute', 1e-12);
