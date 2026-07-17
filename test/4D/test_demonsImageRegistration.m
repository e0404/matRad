function test_suite = test_demonsImageRegistration

test_functions = localfunctions();

initTestSuite;

function test_constructorDefaultsAndStructSerialization
reg = matRad_DemonsImageRegistration();

assertEqual(reg.refScen, 1);
assertEqual(reg.dvfType, 'pull');
assertEqual(reg.dvfUnits, 'mm');
assertEqual(reg.numIterations, 100);
assertEqual(reg.pyramidLevels, 1);
assertEqual(reg.accumulatedFieldSmoothing, 1);

s = struct(reg);
assertEqual(s.className, 'matRad_DemonsImageRegistration');
assertFalse(isfield(s, 'ct'));
assertFalse(isfield(s, 'cst'));

s.dvfType = 'push';
s.numIterations = 5;
restoredReg = matRad_DemonsImageRegistration(s);
assertEqual(restoredReg.dvfType, 'push');
assertEqual(restoredReg.numIterations, 5);
assertEqual(restoredReg.refScen, reg.refScen);
assertEqual(restoredReg.dvfUnits, reg.dvfUnits);

% unknown fields are ignored with a warning, not an error
s.unknownProperty = 42;
restoredReg = matRad_DemonsImageRegistration(s);
assertEqual(restoredReg.numIterations, 5);

function test_invalidPropertyValuesThrow
reg = matRad_DemonsImageRegistration();

assertExceptionThrown(@() helper_assignProperty(reg, 'dvfType', 'invalid'), 'matRad:Error');
assertExceptionThrown(@() helper_assignProperty(reg, 'dvfUnits', 'cm'), 'matRad:Error');
assertExceptionThrown(@() helper_assignProperty(reg, 'refScen', 0), 'matRad:Error');
assertExceptionThrown(@() helper_assignProperty(reg, 'numIterations', 0.5), 'matRad:Error');
assertExceptionThrown(@() helper_assignProperty(reg, 'pyramidLevels', -1), 'matRad:Error');
assertExceptionThrown(@() helper_assignProperty(reg, 'accumulatedFieldSmoothing', -1), 'matRad:Error');

% invalid values also throw when restoring from struct
assertExceptionThrown(@() matRad_DemonsImageRegistration(struct('dvfType', 'invalid')), 'matRad:Error');

function test_calcDvfValidatesCt
[ct, ~] = helper_registrationFixture();

reg = matRad_DemonsImageRegistration();
reg.refScen = 3;
assertExceptionThrown(@() reg.calcDVF(ct), 'matRad:Error');

reg = matRad_DemonsImageRegistration();
ctBroken = ct;
ctBroken.cubeHU{2} = zeros([ct.cubeDim(1) ct.cubeDim(2) ct.cubeDim(3) + 1]);
assertExceptionThrown(@() reg.calcDVF(ctBroken), 'matRad:Error');

ctBroken = rmfield(ct, 'cubeHU');
assertExceptionThrown(@() reg.calcDVF(ctBroken), 'matRad:Error');

function test_propContoursRejectsPullDvf
[ct, cst] = helper_registrationFixture();

% pull-configured registration without existing DVFs
reg = matRad_DemonsImageRegistration(struct('dvfType', 'pull'));
assertExceptionThrown(@() reg.propContours(ct, cst), 'matRad:Error');

% existing pull DVFs traveling with the ct
ct = helper_addZeroDvf(ct, 'pull');
reg = matRad_DemonsImageRegistration(struct('dvfType', 'push'));
assertExceptionThrown(@() reg.propContours(ct, cst), 'matRad:Error');

function test_calcDvfFailsWhenDependenciesMissing
if helper_hasImageRegistrationDependencies()
    moxunit_throw_test_skipped_exception('Image Processing Toolbox is available.');
end

[ct, ~] = helper_registrationFixture();
reg = matRad_DemonsImageRegistration();

assertExceptionThrown(@() reg.calcDVF(ct), 'matRad:Error');

function test_calcDvfZeroMotionPull
if ~helper_hasImageRegistrationDependencies()
    moxunit_throw_test_skipped_exception('Image Processing Toolbox is not available.');
end

[ct, ~] = helper_registrationFixture();
ct.cubeHU{2} = ct.cubeHU{1}; % no motion

reg = matRad_DemonsImageRegistration(struct('numIterations', 10));
ct = reg.calcDVF(ct);

assertEqual(ct.dvfMetadata.dvfType, 'pull');
assertEqual(ct.dvfMetadata.dvfUnits, 'mm');
assertEqual(ct.dvfMetadata.refScen, 1);
assertEqual(size(ct.dvf{2}), [3 ct.cubeDim]);
assertEqual(ct.dvf{1}, zeros([3 ct.cubeDim]));
assertTrue(max(abs(ct.dvf{2}(:))) < 1e-6);

function test_calcDvfPullConventionMatchesDoseAcc
% Registers a gaussian blob shifted by one voxel in x and verifies that the
% resulting pull DVFs accumulate dose correctly via matRad_doseAcc (DDM),
% i.e. sign, component order, and mm units all match matRad's conventions.
if ~helper_hasImageRegistrationDependencies()
    moxunit_throw_test_skipped_exception('Image Processing Toolbox is not available.');
end

[ct, cst] = helper_registrationFixture();

reg = matRad_DemonsImageRegistration();
ct = reg.calcDVF(ct);

% the blob in scenario 2 sits one voxel further in x; the pull field on the
% reference grid must therefore point by about -1 voxel * resolution.x mm
blobMask = ct.cubeHU{1} > 30;
dvfX = squeeze(ct.dvf{2}(1, :, :, :));
assertTrue(mean(dvfX(blobMask)) < -0.5 * ct.resolution.x);

% dose moving with the anatomy must accumulate to twice the reference dose
phaseCubes = {ct.cubeHU{1}, ct.cubeHU{2}};
dAcc = matRad_doseAcc(ct, phaseCubes, cst, 'DDM');

ix = cst{1, 4}{1, 1};
dExpected = 2 * ct.cubeHU{1};
errAcc = norm(dAcc(ix) - dExpected(ix));
errNaive = norm(ct.cubeHU{2}(ix) - ct.cubeHU{1}(ix)); % accumulation without any warping
assertTrue(errAcc < 0.5 * errNaive);

function test_propContoursConstantPushShift
% Constant synthetic push DVF of one voxel in x: since push fields satisfy
% scenCube(pos) = refCube(pos + dvf(pos)), a field of -1 voxel moves the
% contour by +1 voxel. Run in mm (anisotropic resolution) and voxel units
% to also cover the unit conversion.
if ~helper_hasImageRegistrationDependencies()
    moxunit_throw_test_skipped_exception('Image Processing Toolbox is not available.');
end

unitSet = {'mm', 'voxel'};

for u = 1:numel(unitSet)
    [ct, cst] = helper_registrationFixture();
    ct = helper_addZeroDvf(ct, 'push');
    ct.dvfMetadata.dvfUnits = unitSet{u};

    shiftVoxel = -1;
    if strcmp(unitSet{u}, 'mm')
        ct.dvf{2}(1, :, :, :) = shiftVoxel * ct.resolution.x;
    else
        ct.dvf{2}(1, :, :, :) = shiftVoxel;
    end

    reg = matRad_DemonsImageRegistration(struct('dvfType', 'push'));
    [~, cstOut] = reg.propContours(ct, cst);

    [y, x, z] = ind2sub(ct.cubeDim, cst{1, 4}{1, 1});
    expectedContour = sort(sub2ind(ct.cubeDim, y, x + 1, z));

    assertEqual(sort(cstOut{1, 4}{1, 2}), expectedContour);
    assertEqual(sort(cstOut{1, 4}{1, 1}), sort(cst{1, 4}{1, 1})); % reference untouched
end

function test_propContoursPreservesEmptyStructures
if ~helper_hasImageRegistrationDependencies()
    moxunit_throw_test_skipped_exception('Image Processing Toolbox is not available.');
end

[ct, cst] = helper_registrationFixture();
cst{1, 4}{1, 1} = [];
ct = helper_addZeroDvf(ct, 'push');

reg = matRad_DemonsImageRegistration(struct('dvfType', 'push'));
[~, cstOut] = reg.propContours(ct, cst);

assertTrue(isempty(cstOut{1, 4}{1, 1}));
assertTrue(isempty(cstOut{1, 4}{1, 2}));

function test_propContoursPushEndToEnd
% Full chain: propContours without precomputed DVFs triggers a push demons
% registration; the propagated contour must follow the blob by one voxel in x
if ~helper_hasImageRegistrationDependencies()
    moxunit_throw_test_skipped_exception('Image Processing Toolbox is not available.');
end

[ct, cst] = helper_registrationFixture();

reg = matRad_DemonsImageRegistration(struct('dvfType', 'push'));
[ct, cstOut] = reg.propContours(ct, cst);

assertEqual(ct.dvfMetadata.dvfType, 'push');
assertFalse(isempty(cstOut{1, 4}{1, 2}));

[yRef, xRef, zRef] = ind2sub(ct.cubeDim, cstOut{1, 4}{1, 1});
[yProp, xProp, zProp] = ind2sub(ct.cubeDim, cstOut{1, 4}{1, 2});
centroidShift = [mean(yProp) mean(xProp) mean(zProp)] - [mean(yRef) mean(xRef) mean(zRef)];

assertElementsAlmostEqual(centroidShift, [0 1 0], 'absolute', 0.75);

% helper functions

function helper_assignProperty(obj, propName, value)
obj.(propName) = value;

function [ct, cst] = helper_registrationFixture()
% Two-scenario ct with a smooth gaussian blob, shifted by one voxel in x in
% the second scenario. Unequal cube dimensions and anisotropic resolution
% help to catch coordinate-swap and unit-scaling bugs.
ct = struct();
ct.numOfCtScen = 2;
ct.cubeDim = [12 14 8];
ct.resolution = struct('x', 2, 'y', 1, 'z', 3);
ct.cubeHU = {helper_gaussianBlob(ct.cubeDim, [6 7 4]), ...
             helper_gaussianBlob(ct.cubeDim, [6 8 4])};

cst = cell(1, 6);
cst{1, 1} = 1;
cst{1, 2} = 'TEST';
cst{1, 3} = 'TARGET';
cst{1, 4} = cell(1, ct.numOfCtScen);
cst{1, 4}{1, 1} = find(ct.cubeHU{1} > 30);
cst{1, 5} = struct('Visible', true);
cst{1, 6} = {};

function cube = helper_gaussianBlob(cubeDim, center)
[y, x, z] = ndgrid(1:cubeDim(1), 1:cubeDim(2), 1:cubeDim(3));
sigma = 1.5;
cube = 100 * exp(-((y - center(1)).^2 + (x - center(2)).^2 + (z - center(3)).^2) / (2 * sigma^2));

function ct = helper_addZeroDvf(ct, dvfType)
ct.dvfMetadata = struct('dvfType', dvfType, 'dvfUnits', 'mm', 'refScen', 1);
ct.dvf = cell(1, ct.numOfCtScen);
for scen = 1:ct.numOfCtScen
    ct.dvf{scen} = zeros([3 ct.cubeDim]);
end

function tf = helper_hasImageRegistrationDependencies()
matRad_cfg = MatRad_Config.instance();
tf = ~matRad_cfg.isOctave && ...
     matRad_checkEnvImageProcessingRequirements() && ...
     exist('imregdemons', 'file') == 2 && ...
     exist('imwarp', 'file') == 2;
