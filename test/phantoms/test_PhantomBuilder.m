function test_suite = test_PhantomBuilder
% The output should always be test_suite, and the function name the same as
% your file name

test_functions = localfunctions();
initTestSuite;

%% Constructor and box VOIs

function test_constructorAcceptsStructAndVectorResolution
% the ct dimensions are swapped to [y x z], the resolution is not
res = struct('x', 2, 'y', 3, 'z', 5);
builderStruct = matRad_PhantomBuilder([20 30 40], res, 1);
builderVector = matRad_PhantomBuilder([20 30 40], [2 3 5], 1);
builderStruct.addSphericalOAR('S', 3);
builderVector.addSphericalOAR('S', 3);
[ctStruct, ~] = builderStruct.getctcst();
[ctVector, ~] = builderVector.getctcst();

assertEqual(ctStruct.resolution, res);
assertEqual(ctVector.resolution, res);
assertEqual(ctStruct.cubeDim, [30 20 40]);
assertEqual(ctStruct.numOfCtScen, 1);

function test_addBoxTargetAndBoxOAR
builder = helper_createBuilder();
builder.addBoxTarget('BoxTarget', [10 10 10]);
builder.addBoxOAR('BoxOAR', [6 6 6]);
[~, cst] = builder.getctcst();

assertEqual(size(cst, 1), 2);
assertEqual(cst{1, 2}, 'BoxTarget');
assertEqual(cst{1, 3}, 'TARGET');
assertEqual(cst{2, 2}, 'BoxOAR');
assertEqual(cst{2, 3}, 'OAR');
assertFalse(isempty(cst{1, 4}{1}));
assertFalse(isempty(cst{2, 4}{1}));

%% Spherical shell (ring) VOIs

function test_addSphericalShellOARAddsCstEntry
builder = helper_createBuilder();
builder.addSphericalTarget('Target', 5);
builder.addSphericalShellOAR('Ring', 5, 8);
[~, cst] = builder.getctcst();

assertEqual(size(cst, 1), 2);
assertEqual(cst{2, 2}, 'Ring');
assertEqual(cst{2, 3}, 'OAR');
assertEqual(numel(builder.volumes), 2);
assertEqual(builder.volumes{2}.innerRadius, 5);
assertEqual(builder.volumes{2}.radius, 8);

function test_addSphericalShellOARYieldsOuterMinusInnerSphere
% the shell must be exactly the outer sphere minus the inner sphere and must
% not overlap a target of the same radius as the inner radius
builder = helper_createBuilder();
builder.addSphericalTarget('Target', 5);
builder.addSphericalShellOAR('Ring', 5, 8);
[~, cst] = builder.getctcst();

reference = helper_createBuilder();
reference.addSphericalOAR('Big', 8);
reference.addSphericalOAR('Small', 5);
[~, cstRef] = reference.getctcst();

assertEqual(sort(cst{2, 4}{1}), sort(setdiff(cstRef{1, 4}{1}, cstRef{2, 4}{1})));
assertTrue(isempty(intersect(cst{1, 4}{1}, cst{2, 4}{1})));
assertFalse(isempty(cst{2, 4}{1}));

function test_addSphericalShellOARAssignsHU
% getctcst indexes obj.volumes by cst row, so a shell must keep that mapping
% intact and receive its own HU
builder = helper_createBuilder();
builder.addSphericalTarget('Target', 5, 'HU', 50);
builder.addSphericalShellOAR('Ring', 5, 8, 'HU', 100);
[ct, cst] = builder.getctcst();

assertTrue(all(ct.cubeHU{1}(cst{2, 4}{1}) == 100));
assertTrue(all(ct.cubeHU{1}(cst{1, 4}{1}) == 50));

function test_addSphericalShellOARAcceptsObjectivesAndOffset
builder = helper_createBuilder();
builder.addSphericalShellOAR('Ring', 3, 6, ...
                             'objectives', DoseObjectives.matRad_SquaredOverdosing(100, 30), ...
                             'offset', [1 0 0]);
[~, cst] = builder.getctcst();

assertEqual(cst{1, 2}, 'Ring');
assertTrue(isa(cst{1, 6}{1}, 'DoseObjectives.matRad_SquaredOverdosing'));
assertEqual(builder.volumes{1}.offset, [1 0 0]);

function test_addSphericalShellOARRejectsInnerRadiusNotBelowOuter
builder = helper_createBuilder();
assertExceptionThrown(@() builder.addSphericalShellOAR('Ring', 8, 5));
assertExceptionThrown(@() builder.addSphericalShellOAR('Ring', 5, 5));
% a failed call must not leave a volume behind
assertEqual(numel(builder.volumes), 0);

function test_addSphericalShellOARWithZeroInnerRadiusEqualsSphere
builder = helper_createBuilder();
builder.addSphericalShellOAR('Shell', 0, 5);
[~, cstShell] = builder.getctcst();

reference = helper_createBuilder();
reference.addSphericalOAR('Sphere', 5);
[~, cstSphere] = reference.getctcst();

assertEqual(sort(cstShell{1, 4}{1}), sort(cstSphere{1, 4}{1}));

%% Helpers

function builder = helper_createBuilder
builder = matRad_PhantomBuilder([30 30 30], [2 2 2], 1);
