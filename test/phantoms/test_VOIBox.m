function test_suite = test_VOIBox
% The output should always be test_suite, and the function name the same as
% your file name

test_functions = localfunctions();
initTestSuite;

%% Constructor Tests

function test_constructorDefaults
box = matRad_PhantomVOIBox('MyBox', 'OAR', [4 6 8]);
assertEqual(box.name, 'MyBox');
assertEqual(box.type, 'OAR');
assertEqual(box.boxDimensions, [4 6 8]);
assertEqual(box.coordType, 'voxel');
assertEqual(box.offset, [0 0 0]);
assertEqual(box.HU, 0);

function test_constructorCustomParams
box = matRad_PhantomVOIBox('MyBox', 'OAR', [4 6 8], ...
                           'coordType', 'mm', 'offset', [1 2 3], 'HU', 200);
assertEqual(box.coordType, 'mm');
assertEqual(box.offset, [1 2 3]);
assertEqual(box.HU, 200);

function test_constructorInvalidInputs
% Non-positive dimensions
assertExceptionThrown(@() matRad_PhantomVOIBox('MyBox', 'OAR', [-1 2 3]));
assertExceptionThrown(@() matRad_PhantomVOIBox('MyBox', 'OAR', [0 2 3]));
% Wrong number of elements
assertExceptionThrown(@() matRad_PhantomVOIBox('MyBox', 'OAR', [4 6]));
assertExceptionThrown(@() matRad_PhantomVOIBox('MyBox', 'OAR', [4 6 8 10]));
% Non-numeric
assertExceptionThrown(@() matRad_PhantomVOIBox('MyBox', 'OAR', 'big'));
% Invalid coordType
assertExceptionThrown(@() matRad_PhantomVOIBox('MyBox', 'OAR', [4 6 8], 'coordType', 'invalid'));

%% set-method Validation

function test_setMethodsValidation
box = matRad_PhantomVOIBox('MyBox', 'OAR', [4 6 8]);
% coordType: valid update and rejection
box.coordType = 'mm';
assertEqual(box.coordType, 'mm');
assertExceptionThrown(@() helper_assignmentTest(box, 'coordType', 'cube'));
% boxDimensions: rejection of bad values
assertExceptionThrown(@() helper_assignmentTest(box, 'boxDimensions', [-1 2 3]));
assertExceptionThrown(@() helper_assignmentTest(box, 'boxDimensions', [0 2 3]));
assertExceptionThrown(@() helper_assignmentTest(box, 'boxDimensions', [4 6]));

%% initializeParameters - CST structure

function test_initializeParametersCstStructure
ct = helper_createTestCt();
cst = {};
box = matRad_PhantomVOIBox('TestBox', 'OAR', [4 6 8]);
cst = box.initializeParameters(ct, cst);
assertEqual(size(cst, 1), 1);
assertEqual(cst{1, 2}, 'TestBox');
assertEqual(cst{1, 3}, 'OAR');

%% initializeParameters - voxel mode

function test_initializeParametersVoxelCoversAll
% A box larger than the CT must be clipped to all voxels.
ct = helper_createTestCt();
cst = {};
box = matRad_PhantomVOIBox('TestBox', 'OAR', [1000 1000 1000]);
cst = box.initializeParameters(ct, cst);
assertEqual(numel(cst{1, 4}{1}), prod(ct.cubeDim));

function test_initializeParametersVoxelAllIndicesValid
ct = helper_createTestCt();
cst = {};
box = matRad_PhantomVOIBox('TestBox', 'OAR', [4 6 2]);
cst = box.initializeParameters(ct, cst);
idx = cst{1, 4}{1};
assertTrue(numel(idx) > 0);
assertTrue(all(idx >= 1) && all(idx <= prod(ct.cubeDim)));

%% initializeParameters - voxel mode axis-permutation tests
%
% Strategy: odd non-square cubeDim=[9 11 7] puts the geometric centre
% exactly at (row=5, col=6, slice=4).  boxDimensions is in [i j k] order:
%   dims(1) = x / col extent,  dims(2) = y / row extent,  dims(3) = z extent.
%
% dims=[3 1 1]: 3 wide in x (cols 5-7), 1 in y (row 5), 1 in z (slice 4).
% dims=[1 3 1]: 1 in x  (col 6),  3 in y (rows 4-6),    1 in z (slice 4).
% The two resulting voxel sets are completely disjoint, so any x/y swap fails.

function test_initializeParametersVoxelAxisPermutation
cubeDim = [9 11 7];
ct = helper_createTestCt(cubeDim, 1);
cRow = 5;
cCol = 6;
cSlice = 4;   % (cubeDim+1)/2

% --- 1x1x1 box: exactly the centre voxel ----------------------------
cst = {};
box = matRad_PhantomVOIBox('B', 'OAR', [1 1 1]);
cst = box.initializeParameters(ct, cst);
assertEqual(cst{1, 4}{1}, sub2ind(cubeDim, cRow, cCol, cSlice));

% --- [3 1 1]: 3 voxels along x / col --------------------------------
% Expected: row=5, cols 5-7, slice=4
cst = {};
box = matRad_PhantomVOIBox('B', 'OAR', [3 1 1]);
cst = box.initializeParameters(ct, cst);
expected = [sub2ind(cubeDim, cRow, cCol - 1, cSlice), ...
            sub2ind(cubeDim, cRow, cCol,     cSlice), ...
            sub2ind(cubeDim, cRow, cCol + 1, cSlice)];
assertEqual(cst{1, 4}{1}, sort(expected)');

% --- [1 3 1]: 3 voxels along y / row --------------------------------
% Expected: rows 4-6, col=6, slice=4
cst = {};
box = matRad_PhantomVOIBox('B', 'OAR', [1 3 1]);
cst = box.initializeParameters(ct, cst);
expected = [sub2ind(cubeDim, cRow - 1, cCol, cSlice), ...
            sub2ind(cubeDim, cRow,     cCol, cSlice), ...
            sub2ind(cubeDim, cRow + 1, cCol, cSlice)];
assertEqual(cst{1, 4}{1}, sort(expected)');

function test_initializeParametersVoxelOffsetAxis
% offset=[1 0 0] is in [i j k]: x/col direction -> col advances by 1.
% offset=[0 1 0] is y/row direction -> row advances by 1.
% A 1x1x1 box isolates a single voxel so the shift is unambiguous.
cubeDim = [9 11 7];
ct = helper_createTestCt(cubeDim, 1);
cRow = 5;
cCol = 6;
cSlice = 4;

% Offset in x (col)
cst = {};
box = matRad_PhantomVOIBox('B', 'OAR', [1 1 1], 'offset', [1 0 0]);
cst = box.initializeParameters(ct, cst);
assertEqual(cst{1, 4}{1}, sub2ind(cubeDim, cRow, cCol + 1, cSlice));

% Offset in y (row)
cst = {};
box = matRad_PhantomVOIBox('B', 'OAR', [1 1 1], 'offset', [0 1 0]);
cst = box.initializeParameters(ct, cst);
assertEqual(cst{1, 4}{1}, sub2ind(cubeDim, cRow + 1, cCol, cSlice));

%% initializeParameters - mm mode axis-permutation tests
%
% Same cubeDim=[9 11 7] with anisotropic resolution res={x=2,y=3,z=5}.
% The geometric centre maps to world coords (ct.x(6), ct.y(5), ct.z(4)).
%
% dims=[4 0.1 0.1]: extends +-2 mm in x (one full voxel each side) -> cols 5-7.
% dims=[0.1 6 0.1]: extends +-3 mm in y (one full voxel each side) -> rows 4-6.
% Any axis swap produces the opposite column/row expansion and fails.

function test_initializeParametersMmAxisPermutation
cubeDim = [9 11 7];
res = struct('x', 2, 'y', 3, 'z', 5);
ct = helper_createTestCt(cubeDim, res);
cRow = 5;
cCol = 6;
cSlice = 4;

% Tiny box: only the centre voxel
cst = {};
box = matRad_PhantomVOIBox('B', 'OAR', [0.5 0.5 0.5], 'coordType', 'mm');
cst = box.initializeParameters(ct, cst);
assertEqual(cst{1, 4}{1}, sub2ind(cubeDim, cRow, cCol, cSlice));

% dims=[2*res.x, narrow, narrow]: 3 columns, 1 row, 1 slice
cst = {};
box = matRad_PhantomVOIBox('B', 'OAR', [2 * res.x 0.1 0.1], 'coordType', 'mm');
cst = box.initializeParameters(ct, cst);
expected = [sub2ind(cubeDim, cRow, cCol - 1, cSlice), ...
            sub2ind(cubeDim, cRow, cCol,     cSlice), ...
            sub2ind(cubeDim, cRow, cCol + 1, cSlice)];
assertEqual(cst{1, 4}{1}, sort(expected)');

% dims=[narrow, 2*res.y, narrow]: 1 column, 3 rows, 1 slice
cst = {};
box = matRad_PhantomVOIBox('B', 'OAR', [0.1 2 * res.y 0.1], 'coordType', 'mm');
cst = box.initializeParameters(ct, cst);
expected = [sub2ind(cubeDim, cRow - 1, cCol, cSlice), ...
            sub2ind(cubeDim, cRow,     cCol, cSlice), ...
            sub2ind(cubeDim, cRow + 1, cCol, cSlice)];
assertEqual(cst{1, 4}{1}, sort(expected)');

function test_initializeParametersMmOffsetAxis
% A step of res.x in the [i j k] x-direction must advance the column;
% a step of res.y in the y-direction must advance the row.
cubeDim = [9 11 7];
res = struct('x', 2, 'y', 3, 'z', 5);
ct = helper_createTestCt(cubeDim, res);
cRow = 5;
cCol = 6;
cSlice = 4;

% Offset in x (col)
cst = {};
box = matRad_PhantomVOIBox('B', 'OAR', [0.5 0.5 0.5], 'coordType', 'mm', ...
                           'offset', [res.x 0 0]);
cst = box.initializeParameters(ct, cst);
assertEqual(cst{1, 4}{1}, sub2ind(cubeDim, cRow, cCol + 1, cSlice));

% Offset in y (row)
cst = {};
box = matRad_PhantomVOIBox('B', 'OAR', [0.5 0.5 0.5], 'coordType', 'mm', ...
                           'offset', [0 res.y 0]);
cst = box.initializeParameters(ct, cst);
assertEqual(cst{1, 4}{1}, sub2ind(cubeDim, cRow + 1, cCol, cSlice));
