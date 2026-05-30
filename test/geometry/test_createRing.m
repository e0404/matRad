function test_suite = test_createRing

test_functions = localfunctions();

initTestSuite;

function test_createRingBuildsOuterMinusInnerVoi
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin, metadata] = helper_ringArguments();

[cst, ixRing] = matRad_createRing(1, 2, cst, ct, outerMargin, innerMargin, metadata);

targetVoxel = sub2ind(ct.cubeDim, 3, 3, 2);
expectedRing = sort([ ...
                     sub2ind(ct.cubeDim, 2, 3, 2); ...
                     sub2ind(ct.cubeDim, 4, 3, 2); ...
                     sub2ind(ct.cubeDim, 3, 2, 2); ...
                     sub2ind(ct.cubeDim, 3, 4, 2); ...
                     sub2ind(ct.cubeDim, 3, 3, 1); ...
                     sub2ind(ct.cubeDim, 3, 3, 3)]);

assertEqual(ixRing, 3);
assertEqual(cst{ixRing, 2}, 'PTV_RING');
assertEqual(cst{ixRing, 3}, 'OAR');
assertEqual(cst{ixRing, 5}.visibleColor, [0 1 0]);
assertFalse(any(cst{ixRing, 4}{1} == targetVoxel));
assertEqual(sort(cst{ixRing, 4}{1}), expectedRing);

function test_createRingClipsToLimitVoi
[ct, cst] = helper_createRingFixture();
cst{2, 4}{1} = sub2ind(ct.cubeDim, [2 3 4], [3 3 3], [2 2 2])';
[outerMargin, innerMargin, metadata] = helper_ringArguments();

[cst, ixRing] = matRad_createRing(1, 2, cst, ct, outerMargin, innerMargin, metadata);

expectedRing = sort([ ...
                     sub2ind(ct.cubeDim, 2, 3, 2); ...
                     sub2ind(ct.cubeDim, 4, 3, 2)]);

assertEqual(sort(cst{ixRing, 4}{1}), expectedRing);

function [ct, cst] = helper_createRingFixture
ct.cubeDim = [5 5 3];
ct.numOfCtScen = 1;
ct.resolution = struct('x', 1, 'y', 1, 'z', 1);

targetVoxel = sub2ind(ct.cubeDim, 3, 3, 2);
bodyVoxels = (1:prod(ct.cubeDim))';

cst = cell(2, 6);
cst{1, 1} = 0;
cst{1, 2} = 'PTV';
cst{1, 3} = 'TARGET';
cst{1, 4} = {targetVoxel};
cst{1, 5} = struct('visibleColor', [1 0 0]);
cst{1, 6} = {};

cst{2, 1} = 1;
cst{2, 2} = 'BODY';
cst{2, 3} = 'OAR';
cst{2, 4} = {bodyVoxels};
cst{2, 5} = struct('visibleColor', [0 0 1]);
cst{2, 6} = {};

function [outerMargin, innerMargin, metadata] = helper_ringArguments
metadata.name = 'PTV_RING';
metadata.type = 'OAR';
metadata.visibleColor = [0 1 0];
outerMargin = struct('x', 1, 'y', 1, 'z', 1);
innerMargin = struct('x', 0, 'y', 0, 'z', 0);
