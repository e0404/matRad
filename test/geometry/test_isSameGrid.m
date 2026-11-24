function test_suite = test_isSameGrid
% Tests for matRad_isSameGrid, which compares voxel positions and not only
% the dimensions of two grids (#937)

test_functions = localfunctions();

initTestSuite;

function grid = helper_createGrid(res, numOfVoxels)
grid.resolution = struct('x', res(1), 'y', res(2), 'z', res(3));
grid.x = res(1) * (1:numOfVoxels(1));
grid.y = res(2) * (1:numOfVoxels(2));
grid.z = res(3) * (1:numOfVoxels(3));
grid.dimensions = [numOfVoxels(2) numOfVoxels(1) numOfVoxels(3)];
grid.numOfVoxels = prod(grid.dimensions);

function test_identicalGridsAreSame
grid = helper_createGrid([3 3 2.5], [22 10 8]);
assertTrue(matRad_isSameGrid(grid, grid));

function test_differentDimensionsAreNotSame
gridA = helper_createGrid([3 3 3], [22 10 8]);
gridB = helper_createGrid([3 3 3], [21 10 8]);
assertFalse(matRad_isSameGrid(gridA, gridB));

function test_sameDimensionsDifferentAxesAreNotSame
% 22 CT voxels at 3 mm and a centered dose grid at 3.1 mm have the same
% dimensions but different voxel positions
ctGrid = helper_createGrid([3 3 3], [22 10 8]);
doseGrid = ctGrid;
[doseGrid.x, doseGrid.y, doseGrid.z] = matRad_createDoseGridAxes(ctGrid, struct('x', 3.1, 'y', 3, 'z', 3));
doseGrid.resolution.x = 3.1;
assertEqual(numel(doseGrid.x), numel(ctGrid.x));
assertFalse(matRad_isSameGrid(ctGrid, doseGrid));

function test_equalResolutionDoseGridIsSame
ctGrid = helper_createGrid([3 3 2.5], [22 10 8]);
doseGrid = ctGrid;
[doseGrid.x, doseGrid.y, doseGrid.z] = matRad_createDoseGridAxes(ctGrid, ctGrid.resolution);
assertTrue(matRad_isSameGrid(ctGrid, doseGrid));

function test_rowAndColumnAxesAreSame
gridA = helper_createGrid([3 3 3], [22 10 8]);
gridB = gridA;
gridB.x = gridB.x';
assertTrue(matRad_isSameGrid(gridA, gridB));

function test_resolutionIsComparedWithoutAxes
gridA = rmfield(helper_createGrid([3 3 3], [22 10 8]), {'x', 'y', 'z'});
gridB = gridA;
assertTrue(matRad_isSameGrid(gridA, gridB));
gridB.resolution.x = 3.1;
assertFalse(matRad_isSameGrid(gridA, gridB));
