function test_suite = test_elasticImageRegistration

test_functions = localfunctions();

initTestSuite;

function test_constructorDefaultsAndStructSerialization
[ct, cst] = helper_registrationFixture();

registration = matRad_ElasticImageRegistration(ct, cst);
registrationStruct = struct(registration);
restoredRegistration = matRad_ElasticImageRegistration(registrationStruct);

assertEqual(registration.refScen, 1);
assertEqual(registration.metadata.dvfType, 'pull');
assertEqual(registration.metadata.dvfUnits, 'voxel');
assertEqual(registration.metadata.numIterations, 100);
assertEqual(registration.metadata.pyramidLevels, 1);
assertEqual(registration.metadata.accumulatedFieldSmoothing, 1);
assertEqual(registrationStruct.className, 'matRad_ElasticImageRegistration');
assertEqual(restoredRegistration.refScen, registration.refScen);
assertEqual(restoredRegistration.metadata, registration.metadata);

function test_invalidReferenceScenarioThrowsError
[ct, cst] = helper_registrationFixture();

assertExceptionThrown(@() matRad_ElasticImageRegistration(ct, cst, 3), 'matRad:Error');

function test_invalidDvfTypeThrowsError
[ct, cst] = helper_registrationFixture();

assertExceptionThrown(@() matRad_ElasticImageRegistration(ct, cst, 1, struct('dvfType', 'invalid')), ...
                      'matRad:Error');

function test_invalidCtScenarioCubeSizeThrowsError
[ct, cst] = helper_registrationFixture();
ct.cubeHU{2} = zeros([ct.cubeDim(1) ct.cubeDim(2) ct.cubeDim(3) + 1]);

assertExceptionThrown(@() matRad_ElasticImageRegistration(ct, cst), 'matRad:Error');

function test_calcDvfFailsWhenImageProcessingToolboxIsMissing
if helper_hasElasticRegistrationDependencies()
    moxunit_throw_test_skipped_exception('Image Processing Toolbox is available.');
end

[ct, cst] = helper_registrationFixture();
registration = matRad_ElasticImageRegistration(ct, cst, 1, struct('dvfType', 'pull'));

assertExceptionThrown(@() registration.calcDVF(), 'matRad:Error');

function test_calcDvfBuildsPullDvfMetadata
if ~helper_hasElasticRegistrationDependencies()
    moxunit_throw_test_skipped_exception('Image Processing Toolbox is not available.');
end

[ct, cst] = helper_registrationFixture();
metadata = struct('dvfType', 'pull', 'numIterations', 1);
registration = matRad_ElasticImageRegistration(ct, cst, 1, metadata);

ctOut = registration.calcDVF();

assertEqual(ctOut.refScen, 1);
assertEqual(ctOut.dvfType, 'pull');
assertEqual(ctOut.dvfUnits, 'voxel');
assertEqual(ctOut.dvfDim, ct.cubeDim);
assertEqual(ctOut.dvfMetadata.refScen, 1);
assertEqual(ctOut.dvfMetadata.referenceCtScen, 1);
assertEqual(size(ctOut.dvf{1}), [3 ct.cubeDim]);
assertElementsAlmostEqual(ctOut.dvf{1}, zeros([3 ct.cubeDim]), 'absolute', 1e-12);
assertTrue(max(abs(ctOut.dvf{2}(:))) < 1e-4);

function test_propContoursRejectsPullDvf
if ~helper_hasElasticRegistrationDependencies()
    moxunit_throw_test_skipped_exception('Image Processing Toolbox is not available.');
end

[ct, cst] = helper_registrationFixture();
ct = helper_addZeroDvf(ct, 'pull');
registration = matRad_ElasticImageRegistration(ct, cst, 1, struct('dvfType', 'pull'));

assertExceptionThrown(@() registration.propContours(), 'matRad:Error');

function test_propContoursPreservesEmptyStructures
if ~helper_hasElasticRegistrationDependencies()
    moxunit_throw_test_skipped_exception('Image Processing Toolbox is not available.');
end

[ct, cst] = helper_registrationFixture();
cst{1, 4}{1, 1} = [];
ct = helper_addZeroDvf(ct, 'push');
registration = matRad_ElasticImageRegistration(ct, cst, 1, struct('dvfType', 'push'));

[~, cstOut] = registration.propContours();

assertTrue(isempty(cstOut{1, 4}{1, 1}));
assertTrue(isempty(cstOut{1, 4}{1, 2}));

function test_propContoursMapsMaskWithZeroPushDvf
if ~helper_hasElasticRegistrationDependencies()
    moxunit_throw_test_skipped_exception('Image Processing Toolbox is not available.');
end

[ct, cst] = helper_registrationFixture();
ct = helper_addZeroDvf(ct, 'push');
registration = matRad_ElasticImageRegistration(ct, cst, 1, struct('dvfType', 'push'));

[~, cstOut] = registration.propContours();

assertEqual(sort(cstOut{1, 4}{1, 2}(:)), sort(cst{1, 4}{1, 1}(:)));

function [ct, cst] = helper_registrationFixture()
ct = struct();
ct.numOfCtScen = 2;
ct.cubeDim = [5 5 3];
ct.resolution = struct('x', 1, 'y', 1, 'z', 1);

cube = zeros(ct.cubeDim);
cube(3, 3, 2) = 100;
ct.cubeHU = {cube, cube};

cst = cell(1, 6);
cst{1, 1} = 1;
cst{1, 2} = 'TEST';
cst{1, 3} = 'TARGET';
cst{1, 4} = cell(1, ct.numOfCtScen);
cst{1, 4}{1, 1} = sub2ind(ct.cubeDim, [3 3], [3 4], [2 2]);
cst{1, 5} = struct('Visible', true);
cst{1, 6} = {};

function ct = helper_addZeroDvf(ct, dvfType)
ct.refScen = 1;
ct.dvfType = dvfType;
ct.dvfUnits = 'voxel';
ct.dvfDim = ct.cubeDim;
ct.dvfMetadata.dvfType = dvfType;
ct.dvfMetadata.dvfUnits = 'voxel';
ct.dvfMetadata.refScen = 1;
ct.dvfMetadata.referenceCtScen = 1;
ct.dvf = cell(1, ct.numOfCtScen);
for scen = 1:ct.numOfCtScen
    ct.dvf{scen} = zeros([3 ct.cubeDim]);
end

function tf = helper_hasElasticRegistrationDependencies()
env = matRad_getEnvironment();
tf = strcmp(env, 'MATLAB') && exist('imregdemons', 'file') == 2 && exist('imwarp', 'file') == 2;
if tf && (exist('license', 'builtin') == 5 || exist('license', 'file') == 2)
    try
        tf = license('test', 'image_toolbox');
    catch
    end
end
