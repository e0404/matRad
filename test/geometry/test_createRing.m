function test_suite = test_createRing

test_functions = localfunctions();

initTestSuite;

function test_createRingBuildsOuterMinusInnerVoi
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin, metadata] = helper_ringArguments();

[cst, ixRing] = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2);

targetVoxel = sub2ind(ct.cubeDim, 3, 3, 2);
expectedRing = sort([ ...
                     sub2ind(ct.cubeDim, 2, 3, 2); ...
                     sub2ind(ct.cubeDim, 4, 3, 2); ...
                     sub2ind(ct.cubeDim, 3, 2, 2); ...
                     sub2ind(ct.cubeDim, 3, 4, 2); ...
                     sub2ind(ct.cubeDim, 3, 3, 1); ...
                     sub2ind(ct.cubeDim, 3, 3, 3)]);

assertEqual(ixRing, 3);
assertEqual(cst{ixRing, 1}, 2);
assertEqual(cst{ixRing, 2}, 'PTV_RING');
assertEqual(cst{ixRing, 3}, 'OAR');
assertEqual(cst{ixRing, 5}.visibleColor, [0 1 0]);
assertEqual(cst{ixRing, 6}, {});
assertFalse(any(cst{ixRing, 4}{1} == targetVoxel));
assertEqual(sort(cst{ixRing, 4}{1}), expectedRing);

function test_createRingClipsToLimitVoi
[ct, cst] = helper_createRingFixture();
cst{2, 4}{1} = sub2ind(ct.cubeDim, [2 3 4], [3 3 3], [2 2 2])';
[outerMargin, innerMargin, metadata] = helper_ringArguments();

[cst, ixRing] = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2);

expectedRing = sort([ ...
                     sub2ind(ct.cubeDim, 2, 3, 2); ...
                     sub2ind(ct.cubeDim, 4, 3, 2)]);

assertEqual(sort(cst{ixRing, 4}{1}), expectedRing);

function test_createRingDoesNotGrowThroughVoisOutsideTheLimit
% the margin must not travel through a structure outside the limiting VOI
% and re-enter it further out
[ct, cst] = helper_createRingFixture();

corridorVoxel = sub2ind(ct.cubeDim, 3, 4, 2);
leakVoxel     = sub2ind(ct.cubeDim, 3, 5, 2);

cst{2, 4}{1} = setdiff((1:prod(ct.cubeDim))', corridorVoxel);
cst{3, 1} = 2;
cst{3, 2} = 'CORRIDOR';
cst{3, 3} = 'OAR';
cst{3, 4} = {corridorVoxel};
cst{3, 5} = struct('visibleColor', [1 1 0], 'Priority', 3);
cst{3, 6} = {};

[~, innerMargin, metadata] = helper_ringArguments();
outerMargin = struct('x', 2, 'y', 2, 'z', 2);

[cst, ixRing] = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2);

assertFalse(any(cst{ixRing, 4}{1} == corridorVoxel));
assertFalse(any(cst{ixRing, 4}{1} == leakVoxel));

function test_createRingHandlesMultipleCtScenarios
[ct, cst] = helper_createRingFixture();
ct.numOfCtScen = 2;
cst{1, 4} = {sub2ind(ct.cubeDim, 3, 3, 2), sub2ind(ct.cubeDim, 2, 2, 2)};
cst{2, 4} = {(1:prod(ct.cubeDim))', (1:prod(ct.cubeDim))'};
[outerMargin, innerMargin, metadata] = helper_ringArguments();

[cst, ixRing] = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2);

assertEqual(numel(cst{ixRing, 4}), 2);
assertEqual(sort(cst{ixRing, 4}{1}), sort([ ...
                                           sub2ind(ct.cubeDim, 2, 3, 2); ...
                                           sub2ind(ct.cubeDim, 4, 3, 2); ...
                                           sub2ind(ct.cubeDim, 3, 2, 2); ...
                                           sub2ind(ct.cubeDim, 3, 4, 2); ...
                                           sub2ind(ct.cubeDim, 3, 3, 1); ...
                                           sub2ind(ct.cubeDim, 3, 3, 3)]));
assertEqual(sort(cst{ixRing, 4}{2}), sort([ ...
                                           sub2ind(ct.cubeDim, 1, 2, 2); ...
                                           sub2ind(ct.cubeDim, 3, 2, 2); ...
                                           sub2ind(ct.cubeDim, 2, 1, 2); ...
                                           sub2ind(ct.cubeDim, 2, 3, 2); ...
                                           sub2ind(ct.cubeDim, 2, 2, 1); ...
                                           sub2ind(ct.cubeDim, 2, 2, 3)]));

function test_createRingUnionsBaseVoiAcrossScenariosOnRequest
% with bUnionScen the ring is grown around the union of the base VOI, but
% the limiting VOI still applies per scenario
[ct, cst] = helper_createRingFixture();
ct.numOfCtScen = 2;
cst{1, 4} = {sub2ind(ct.cubeDim, 3, 3, 2), sub2ind(ct.cubeDim, 2, 2, 2)};

clippedVoxel = sub2ind(ct.cubeDim, 1, 2, 2);
cst{2, 4} = {(1:prod(ct.cubeDim))', setdiff((1:prod(ct.cubeDim))', clippedVoxel)};

[outerMargin, innerMargin, metadata] = helper_ringArguments();

[cst, ixRing] = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2, 'unionScenarios', true);

% neighbours of both base voxels, minus the two base voxels themselves
expectedRing = sort([ ...
                     sub2ind(ct.cubeDim, 2, 3, 2); ...
                     sub2ind(ct.cubeDim, 4, 3, 2); ...
                     sub2ind(ct.cubeDim, 3, 2, 2); ...
                     sub2ind(ct.cubeDim, 3, 4, 2); ...
                     sub2ind(ct.cubeDim, 3, 3, 1); ...
                     sub2ind(ct.cubeDim, 3, 3, 3); ...
                     sub2ind(ct.cubeDim, 1, 2, 2); ...
                     sub2ind(ct.cubeDim, 2, 1, 2); ...
                     sub2ind(ct.cubeDim, 2, 2, 1); ...
                     sub2ind(ct.cubeDim, 2, 2, 3)]);

assertEqual(sort(cst{ixRing, 4}{1}), expectedRing);
% the second scenario sees the same union base but its own limiting VOI
assertEqual(sort(cst{ixRing, 4}{2}), setdiff(expectedRing, clippedVoxel));

function test_createRingWithoutUnionKeepsScenariosIndependent
% default behaviour must stay per scenario
[ct, cst] = helper_createRingFixture();
ct.numOfCtScen = 2;
cst{1, 4} = {sub2ind(ct.cubeDim, 3, 3, 2), sub2ind(ct.cubeDim, 2, 2, 2)};
cst{2, 4} = {(1:prod(ct.cubeDim))', (1:prod(ct.cubeDim))'};
[outerMargin, innerMargin, metadata] = helper_ringArguments();

[cstUnion, ixRing] = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2, 'unionScenarios', true);
[cstScen, ~] = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2, 'unionScenarios', false);

% per scenario the two rings differ, the union ring is identical everywhere
assertFalse(isequal(sort(cstScen{ixRing, 4}{1}), sort(cstScen{ixRing, 4}{2})));
assertEqual(sort(cstUnion{ixRing, 4}{1}), sort(cstUnion{ixRing, 4}{2}));

function test_createRingRejectsDuplicateVoiName
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin, metadata] = helper_ringArguments();
metadata.name = 'BODY';

assertExceptionThrown(@() matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2));

function test_createRingInheritsPriorityFromBaseVoi
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin, metadata] = helper_ringArguments();

[cst, ixRing] = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2);

assertEqual(cst{ixRing, 5}.Priority, cst{1, 5}.Priority);

function test_createRingUsesPriorityFromMetadata
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin, metadata] = helper_ringArguments();
metadata.Priority = 5;

[cst, ixRing] = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2);

assertEqual(cst{ixRing, 5}.Priority, 5);

function test_createRingStoresObjectivesFromMetadata
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin, metadata] = helper_ringArguments();
metadata.objectives = DoseObjectives.matRad_SquaredOverdosing(100, 30);

[cst, ixRing] = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2);

assertTrue(iscell(cst{ixRing, 6}));
assertEqual(numel(cst{ixRing, 6}), 1);
assertTrue(isa(cst{ixRing, 6}{1}, 'DoseObjectives.matRad_SquaredOverdosing'));

function test_createRingRejectsNegativeMargins
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin, metadata] = helper_ringArguments();
innerMargin.x = -1;

assertExceptionThrown(@() matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2));

function test_createRingRejectsInnerMarginLargerThanOuterMargin
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin, metadata] = helper_ringArguments();
innerMargin.z = 2;

assertExceptionThrown(@() matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2));

function test_createRingRejectsInvalidVoiIndices
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin, metadata] = helper_ringArguments();

assertExceptionThrown(@() matRad_createRing(ct, cst, 0, outerMargin, innerMargin, metadata, 'voiLimit', 2));
assertExceptionThrown(@() matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 7));

function test_createRingAppliesDefaultsForOmittedMetadata
% without any metadata the ring is named after the base VOI and is an OAR
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin] = helper_ringArguments();

[cst, ixRing] = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, 'voiLimit', 2);

assertEqual(cst{ixRing, 2}, 'PTV_RING');
assertEqual(cst{ixRing, 3}, 'OAR');
assertEqual(cst{ixRing, 5}.visibleColor, [0 1 0]);
assertEqual(cst{ixRing, 6}, {});

function test_createRingAcceptsStructAndNameValueInterchangeably
% a metadata struct and the equivalent Name/Value pairs must agree
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin, metadata] = helper_ringArguments();

cstStruct = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2);
cstPairs  = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, 'voiLimit', 2, ...
                              'name', metadata.name, 'type', metadata.type, ...
                              'visibleColor', metadata.visibleColor);

assertEqual(cstStruct(3, :), cstPairs(3, :));

function test_createRingAcceptsStructMixedWithNameValue
% Name/Value pairs after a struct must be applied on top of it
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin, metadata] = helper_ringArguments();

[cst, ixRing] = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2, ...
                                  'type', 'TARGET', 'Priority', 7);

assertEqual(cst{ixRing, 2}, metadata.name);   % from the struct
assertEqual(cst{ixRing, 3}, 'TARGET');        % from the pair
assertEqual(cst{ixRing, 5}.Priority, 7);      % from the pair

function test_createRingAcceptsCstColumnFiveAsMetadata
% the metadata struct mirrors the fifth cst column, so an existing one must
% be usable directly
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin] = helper_ringArguments();

metadata = cst{1, 5};             % visibleColor + Priority
metadata.alphaX = 0.5;
metadata.TissueClass = 3;
metadata.name = 'FROM_CST';

[cst, ixRing] = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2);

assertEqual(cst{ixRing, 2}, 'FROM_CST');
assertEqual(cst{ixRing, 5}.alphaX, 0.5);
assertEqual(cst{ixRing, 5}.TissueClass, 3);
assertEqual(cst{ixRing, 5}.Priority, cst{1, 5}.Priority);

function test_createRingAcceptsVoiNamesAndIndicesInterchangeably
% voiBase and voiLimit accept a cst row index or a VOI name
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin] = helper_ringArguments();

cstByIx   = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, 'voiLimit', 2);
cstByName = matRad_createRing(ct, cst, 'PTV', outerMargin, innerMargin, 'voiLimit', 'BODY');

assertEqual(cstByIx(3, :), cstByName(3, :));

function test_createRingRejectsUnknownVoiName
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin] = helper_ringArguments();

assertExceptionThrown(@() matRad_createRing(ct, cst, 'NOT_A_VOI', outerMargin, innerMargin));
assertExceptionThrown(@() matRad_createRing(ct, cst, 'PTV', outerMargin, innerMargin, 'voiLimit', 'NOT_A_VOI'));

function test_createRingRejectsAmbiguousVoiName
% a duplicate name in the cst cannot be resolved to a single row
[ct, cst] = helper_createRingFixture();
cst{3, 1} = 2;
cst{3, 2} = 'PTV';
cst{3, 3} = 'TARGET';
cst{3, 4} = {sub2ind(ct.cubeDim, 1, 1, 1)};
cst{3, 5} = struct('visibleColor', [1 1 0], 'Priority', 3);
cst{3, 6} = {};
[outerMargin, innerMargin] = helper_ringArguments();

assertExceptionThrown(@() matRad_createRing(ct, cst, 'PTV', outerMargin, innerMargin));

function test_createRingLimitsToAllVoisByDefault
% without voiLimit the ring is clipped to the union of all VOIs, i.e. the
% patient outline. Here BODY only covers a slab, so the ring cannot leave it
[ct, cst] = helper_createRingFixture();
cst{2, 4}{1} = sub2ind(ct.cubeDim, [2 3 4], [3 3 3], [2 2 2])';
[outerMargin, innerMargin] = helper_ringArguments();

[cst, ixRing] = matRad_createRing(ct, cst, 1, outerMargin, innerMargin);

% the union of PTV and BODY is the slab, so only its two free voxels remain
expectedRing = sort([ ...
                     sub2ind(ct.cubeDim, 2, 3, 2); ...
                     sub2ind(ct.cubeDim, 4, 3, 2)]);

assertEqual(sort(cst{ixRing, 4}{1}), expectedRing);

function test_createRingRejectsUnknownParameter
% an unknown option must not be swallowed silently, neither as a Name/Value
% pair nor as a struct field
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin, metadata] = helper_ringArguments();

assertExceptionThrown(@() matRad_createRing(ct, cst, 1, outerMargin, innerMargin, 'voiLimit', 2, 'notAParameter', true));

metadata.notAParameter = true;
assertExceptionThrown(@() matRad_createRing(ct, cst, 1, outerMargin, innerMargin, metadata, 'voiLimit', 2));

function test_createRingRejectsAbbreviatedParameter
% PartialMatching is disabled so that an abbreviation is rejected on Matlab
% and Octave alike, where the inputParser defaults would otherwise differ
[ct, cst] = helper_createRingFixture();
[outerMargin, innerMargin] = helper_ringArguments();

assertExceptionThrown(@() matRad_createRing(ct, cst, 1, outerMargin, innerMargin, 'voiLimit', 2, 'unionScenario', true));
% the full name still works
[~, ixRing] = matRad_createRing(ct, cst, 1, outerMargin, innerMargin, 'voiLimit', 2, 'unionScenarios', true);
assertEqual(ixRing, 3);

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
cst{1, 5} = struct('visibleColor', [1 0 0], 'Priority', 1);
cst{1, 6} = {};

cst{2, 1} = 1;
cst{2, 2} = 'BODY';
cst{2, 3} = 'OAR';
cst{2, 4} = {bodyVoxels};
cst{2, 5} = struct('visibleColor', [0 0 1], 'Priority', 2);
cst{2, 6} = {};

function [outerMargin, innerMargin, metadata] = helper_ringArguments
metadata.name = 'PTV_RING';
metadata.type = 'OAR';
metadata.visibleColor = [0 1 0];
outerMargin = struct('x', 1, 'y', 1, 'z', 1);
innerMargin = struct('x', 0, 'y', 0, 'z', 0);
