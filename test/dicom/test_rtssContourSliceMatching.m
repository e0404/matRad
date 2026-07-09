function test_suite = test_rtssContourSliceMatching

test_functions = localfunctions();

initTestSuite;

function test_acceptsBoundaryContourInsideSliceSlab
ct = helper_ctFixture();

sliceIndices = matRad_findRtssContourSlicesInCt(91.530000000000, ct);

assertEqual(sliceIndices, 3);

function test_rejectsContourOutsideSliceSlab
ct = helper_ctFixture();
outsideZ = ct.z(end) + ct.dicomInfo.SliceThickness / 2 + 0.01;

sliceIndices = matRad_findRtssContourSlicesInCt(outsideZ, ct);

assertTrue(isempty(sliceIndices));

function test_invalidSliceThicknessThrowsError
ct = helper_ctFixture();

ct.dicomInfo.SliceThickness = [];
assertExceptionThrown(@() matRad_findRtssContourSlicesInCt(0, ct), 'matRad:Error');

ct.dicomInfo.SliceThickness = [3 3];
assertExceptionThrown(@() matRad_findRtssContourSlicesInCt(0, ct), 'matRad:Error');

ct.dicomInfo.SliceThickness = '3';
assertExceptionThrown(@() matRad_findRtssContourSlicesInCt(0, ct), 'matRad:Error');

ct.dicomInfo.SliceThickness = 0;
assertExceptionThrown(@() matRad_findRtssContourSlicesInCt(0, ct), 'matRad:Error');

function test_voxelizesBoundaryContourOnExtremeSlice
ct = helper_ctFixture();
structure = helper_structureFixture(91.530000000000);

indices = matRad_convRtssContours2Indices(structure, ct);
[~, ~, iz] = ind2sub(ct.cubeDim, indices);

assertTrue(~isempty(indices));
assertEqual(unique(iz), 3);

function test_voxelizationOmitsContourOutsideCt
ct = helper_ctFixture();
outsideZ = ct.z(end) + ct.dicomInfo.SliceThickness / 2 + 0.01;
structure = helper_structureFixture(outsideZ);

indices = matRad_convRtssContours2Indices(structure, ct);

assertTrue(isempty(indices));

function ct = helper_ctFixture()
ct.x = 1:5;
ct.y = 1:5;
ct.z = [-91.527328693131 0 91.527328693131];
ct.cubeDim = [5 5 3];
ct.cubeHU = {zeros(ct.cubeDim)};
ct.dicomInfo.SlicePositions = ct.z;
ct.dicomInfo.SliceThickness = 3;

function structure = helper_structureFixture(contourZ)
structure.structName = 'BODY';
structure.item(1).points = [2 2 contourZ; ...
                            4 2 contourZ; ...
                            4 4 contourZ; ...
                            2 4 contourZ; ...
                            2 2 contourZ];
