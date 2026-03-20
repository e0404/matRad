function test_suite = test_VOISphere
% The output should always be test_suite, and the function name the same as
% your file name

test_functions = localfunctions();
initTestSuite;

%% Constructor Tests

function test_constructorDefaults
% All default properties from a minimal construction
sphere = matRad_PhantomVOISphere('MySphere', 'OAR', 5);
assertEqual(sphere.name, 'MySphere');
assertEqual(sphere.type, 'OAR');
assertEqual(sphere.radius, 5);
assertEqual(sphere.coordType, 'voxel');
assertEqual(sphere.offset, [0 0 0]);
assertEqual(sphere.HU, 0);

function test_constructorCustomParams
% Custom name-value arguments are stored correctly
sphere = matRad_PhantomVOISphere('MySphere', 'OAR', 5, ...
                                 'coordType', 'mm', 'offset', [1 2 3], 'HU', 200);
assertEqual(sphere.coordType, 'mm');
assertEqual(sphere.offset, [1 2 3]);
assertEqual(sphere.HU, 200);

function test_constructorInvalidInputs
% Invalid radius values
assertExceptionThrown(@() matRad_PhantomVOISphere('MySphere', 'OAR', -1));
assertExceptionThrown(@() matRad_PhantomVOISphere('MySphere', 'OAR', 0));
assertExceptionThrown(@() matRad_PhantomVOISphere('MySphere', 'OAR', [1 2]));
assertExceptionThrown(@() matRad_PhantomVOISphere('MySphere', 'OAR', 'big'));
% Invalid coordType
assertExceptionThrown(@() matRad_PhantomVOISphere('MySphere', 'OAR', 5, 'coordType', 'invalid'));

%% set-method Validation

function test_setMethodsValidation
sphere = matRad_PhantomVOISphere('MySphere', 'OAR', 5);
% coordType: valid update and rejection of invalid value
sphere.coordType = 'mm';
assertEqual(sphere.coordType, 'mm');
assertExceptionThrown(@() helper_assignmentTest(sphere, 'coordType', 'cube'));
% radius: rejection of non-positive and non-scalar values
assertExceptionThrown(@() helper_assignmentTest(sphere, 'radius', -3));
assertExceptionThrown(@() helper_assignmentTest(sphere, 'radius', 0));
assertExceptionThrown(@() helper_assignmentTest(sphere, 'radius', [1 2 3]));

%% initializeParameters Tests

function test_initializeParametersCstStructure
% A single call must add exactly one correctly labelled CST entry
ct = helper_createTestCt();
cst = {};
sphere = matRad_PhantomVOISphere('TestSphere', 'OAR', 3);
cst = sphere.initializeParameters(ct, cst);
assertEqual(size(cst, 1), 1);
assertEqual(cst{1, 2}, 'TestSphere');
assertEqual(cst{1, 3}, 'OAR');

function test_initializeParametersVoxelMode
ct = helper_createTestCt();
cst = {};
radius = 3;
sphere = matRad_PhantomVOISphere('TestSphere', 'OAR', radius);
cst = sphere.initializeParameters(ct, cst);
voxelIndices = cst{1, 4}{1};

% The implementation centers at (cubeDim+1)/2 in [y x z] order.
% Round to the nearest integer voxel for the membership check.
centerPoint = round((ct.cubeDim + 1) / 2);  % [y x z] center voxel
centerLinIx = sub2ind(ct.cubeDim, centerPoint(1), centerPoint(2), centerPoint(3));
assertTrue(ismember(centerLinIx, voxelIndices));

% Corner [1 1 1] is > 7 voxels from the center -> excluded for radius=3
assertFalse(ismember(sub2ind(ct.cubeDim, 1, 1, 1), voxelIndices));

% Voxel count should be close to the continuous sphere volume (4/3)*pi*r^3
expectedVol = (4 / 3) * pi * radius^3;
assertTrue(numel(voxelIndices) > 0.7 * expectedVol);
assertTrue(numel(voxelIndices) < 1.3 * expectedVol);

% All returned indices must be valid linear indices for the cube
assertTrue(all(voxelIndices >= 1) && all(voxelIndices <= prod(ct.cubeDim)));

function test_initializeParametersVoxelLargeRadiusCoversAll
% Center is at (cubeDim+1)/2 = [5.5 6.5 4.5]; the farthest corner is
% sqrt(4.5^2+5.5^2+3.5^2) = ~7.92 voxels away, so radius=9 must capture
% every voxel in the grid.
ct = helper_createTestCt();
cst = {};
sphere = matRad_PhantomVOISphere('TestSphere', 'OAR', 9);
cst = sphere.initializeParameters(ct, cst);
assertEqual(numel(cst{1, 4}{1}), prod(ct.cubeDim));

function test_initializeParametersMmMode
% A larger mm radius must yield strictly more voxels, all within bounds.
% Using non-isotropic resolution means the two radii map to clearly
% different physical volumes even across the coarsest (z=3mm) axis.
ct = helper_createTestCt();
cst = {};
sphereSmall = matRad_PhantomVOISphere('S1', 'OAR',  5, 'coordType', 'mm');
sphereLarge = matRad_PhantomVOISphere('S2', 'OAR', 10, 'coordType', 'mm');
cst1 = sphereSmall.initializeParameters(ct, cst);
cst2 = sphereLarge.initializeParameters(ct, cst);
nSmall = numel(cst1{1, 4}{1});
nLarge = numel(cst2{1, 4}{1});
assertTrue(nSmall > 0);
assertTrue(nLarge > nSmall);
assertTrue(all(cst1{1, 4}{1} >= 1) && all(cst1{1, 4}{1} <= prod(ct.cubeDim)));

function test_initializeParametersLargerRadiusMoreVoxels
% Increasing the radius must strictly grow the voxel set (voxel mode)
ct = helper_createTestCt();
cst = {};
sphereSmall = matRad_PhantomVOISphere('Small', 'OAR', 2);
sphereLarge = matRad_PhantomVOISphere('Large', 'OAR', 4);
cst1 = sphereSmall.initializeParameters(ct, cst);
cst2 = sphereLarge.initializeParameters(ct, cst);
assertTrue(numel(cst2{1, 4}{1}) > numel(cst1{1, 4}{1}));

%% Permutation / axis-direction tests
%
% Strategy: odd non-square cubeDim=[9 11 7] makes (cubeDim+1)/2=[5 6 4]
% land exactly on voxel (row=5, col=6, slice=4).  A radius < 1 voxel
% selects only that single voxel, so an x/y axis swap produces a
% *different linear index* that assertEqual catches immediately.
%
% Offset convention: [i j k] = [x/col, y/row, z/slice].

function test_initializeParametersVoxelAxisPermutation
% In voxel mode an offset of [1 0 0] must advance the column (x/i),
% and [0 1 0] must advance the row (y/j).  A coordinate swap would
% exchange the two, producing the wrong linear index.
cubeDim = [9 11 7];
ct = helper_createTestCt(cubeDim, 1);
cRow = 5;
cCol = 6;
cSlice = 4;   % (cubeDim+1)/2
r = 0.6;      % < 1 voxel: isolates the single centre voxel

% No offset: centre voxel only
cst = {};
cst = matRad_PhantomVOISphere('S', 'OAR', r).initializeParameters(ct, cst);
assertEqual(cst{1, 4}{1}, sub2ind(cubeDim, cRow, cCol, cSlice));

% Offset +1 in i (x / col): col must increase by 1, row unchanged
cst = {};
cst = matRad_PhantomVOISphere('S', 'OAR', r, 'offset', [1 0 0]).initializeParameters(ct, cst);
assertEqual(cst{1, 4}{1}, sub2ind(cubeDim, cRow, cCol + 1, cSlice));

% Offset +1 in j (y / row): row must increase by 1, col unchanged
cst = {};
cst = matRad_PhantomVOISphere('S', 'OAR', r, 'offset', [0 1 0]).initializeParameters(ct, cst);
assertEqual(cst{1, 4}{1}, sub2ind(cubeDim, cRow + 1, cCol, cSlice));

function test_initializeParametersMmAxisPermutation
% Same idea in mm mode.  An offset of exactly one voxel spacing in x
% must advance the column index; one spacing in y must advance the row.
% Anisotropic resolution (x!=y) makes a coordinate swap detectable even
% if voxel counts happen to be equal.
cubeDim = [9 11 7];
res = struct('x', 2, 'y', 3, 'z', 5);
ct = helper_createTestCt(cubeDim, res);
cRow = 5;
cCol = 6;
cSlice = 4;   % (cubeDim+1)/2
r = 0.8;      % < min(res) = 2 mm: isolates centre voxel

% No offset: centre voxel only
cst = {};
cst = matRad_PhantomVOISphere('S', 'OAR', r, 'coordType', 'mm').initializeParameters(ct, cst);
assertEqual(cst{1, 4}{1}, sub2ind(cubeDim, cRow, cCol, cSlice));

% Offset one x-spacing in i (x / col): col must increase by 1
cst = {};
cst = matRad_PhantomVOISphere('S', 'OAR', r, 'coordType', 'mm', ...
                              'offset', [res.x 0 0]).initializeParameters(ct, cst);
assertEqual(cst{1, 4}{1}, sub2ind(cubeDim, cRow, cCol + 1, cSlice));

% Offset one y-spacing in j (y / row): row must increase by 1
cst = {};
cst = matRad_PhantomVOISphere('S', 'OAR', r, 'coordType', 'mm', ...
                              'offset', [0 res.y 0]).initializeParameters(ct, cst);
assertEqual(cst{1, 4}{1}, sub2ind(cubeDim, cRow + 1, cCol, cSlice));
