function test_suite = test_isgpuarrayCompat
% The output should always be test_suite, and the function name the same as
% your file name

% To collect all tests defined below, this is needed in newer Matlab
% versions. test_functions will collect function handles to below test
% functions
test_functions = localfunctions();

% This will initialize the test suite, i.e., take the functions from
% test_functions, check if they contain "test", convert them into a MOxUnit
% Test Case, and add them to the test-runner
initTestSuite;

function test_isgpuarrayCompatRegularArrays
assertFalse(matRad_isgpuarrayCompat(1));
assertFalse(matRad_isgpuarrayCompat([1, 2, 3]));
assertFalse(matRad_isgpuarrayCompat(single(ones(3, 2))));
assertFalse(matRad_isgpuarrayCompat(int32(4)));
assertFalse(matRad_isgpuarrayCompat(true));
assertFalse(matRad_isgpuarrayCompat([]));
assertFalse(matRad_isgpuarrayCompat('text'));
assertFalse(matRad_isgpuarrayCompat({1, 2}));
assertFalse(matRad_isgpuarrayCompat(struct('a', 1)));
assertFalse(matRad_isgpuarrayCompat(sparse([1, 0, 2])));

function test_isgpuarrayCompatReturnsLogicalScalar
tf = matRad_isgpuarrayCompat(ones(4));
assertTrue(islogical(tf));
assertTrue(isscalar(tf));

function test_isgpuarrayCompatDoesNotNeedIsgpuarray
% the wrapper must not rely on isgpuarray, which requires the Parallel
% Computing Toolbox (#938)
cleanupObj = helper_shadowIsgpuarray(); %#ok<NASGU>
assertFalse(matRad_isgpuarrayCompat(ones(3)));

function test_isgpuarrayCompatGpuArray
matRad_cfg = MatRad_Config.instance();
if matRad_cfg.isOctave
    moxunit_throw_test_skipped_exception('gpuArray does not exist in Octave');
end
if ~helper_gpuAvailable()
    moxunit_throw_test_skipped_exception('No Parallel Computing Toolbox or supported GPU available');
end

g = gpuArray(ones(3));
assertTrue(matRad_isgpuarrayCompat(g));
assertFalse(matRad_isgpuarrayCompat(gather(g)));

function test_interp1GriddedInterpolantWithoutIsgpuarray
% matRad_interp1 used isgpuarray on its griddedInterpolant branch, which
% failed without the Parallel Computing Toolbox (#938)
cleanupObj = helper_shadowIsgpuarray(); %#ok<NASGU>

xi = (0:10)';
yi = 2 * xi;
x = [0.5; 3.25; 7.75];
% more than one query point selects the griddedInterpolant branch in Matlab
y = matRad_interp1(xi, yi, x);
assertElementsAlmostEqual(y, 2 * x);

y = matRad_interp1(xi, [yi, 3 * yi], x);
assertElementsAlmostEqual(y, [2 * x, 6 * x]);

function test_rayTracerTraceCubeWithoutIsgpuarray
% matRad_RayTracer.traceCube used isgpuarray, which failed without the
% Parallel Computing Toolbox (#938)
cleanupObj = helper_shadowIsgpuarray(); %#ok<NASGU>

ct = helper_createTestCt([8, 8, 8], 2, 'createCoordinateArrays', true);
cubes = {ones(ct.cubeDim)};
rt = matRad_RayTracerSiddon(cubes, ct);

stfElement.gantryAngle = 0;
stfElement.couchAngle = 0;
stfElement.isoCenter = [0, 0, 0];
stfElement.SAD = 1000;
stfElement.sourcePoint = [0, -stfElement.SAD, 0];
stfElement.sourcePoint_bev = [0, -stfElement.SAD, 0];
stfElement.ray(1).rayPos_bev = [0, 0, 0];

radDepthsV = rt.traceCube(stfElement);
assertTrue(iscell(radDepthsV));
assertEqual(size(radDepthsV{1}), [prod(ct.cubeDim), 1]);
assertTrue(any(isfinite(radDepthsV{1})));

function cleanupObj = helper_shadowIsgpuarray()
% puts an isgpuarray function that always errors in front of the Matlab path,
% mimicking a Matlab installation without the Parallel Computing Toolbox.
% The shadow is removed when the returned onCleanup object is destroyed.
shadowDir = tempname();
mkdir(shadowDir);
fid = fopen(fullfile(shadowDir, 'isgpuarray.m'), 'w');
fprintf(fid, 'function tf = isgpuarray(varargin)\n');
fprintf(fid, 'error(''matRad:test:isgpuarrayShadow'', ''isgpuarray must not be called'');\n');
fprintf(fid, 'end\n');
fclose(fid);
addpath(shadowDir);
% make sure that the shadow takes effect for already cached lookups
rehash;
cleanupObj = onCleanup(@() helper_removeShadow(shadowDir));

function helper_removeShadow(shadowDir)
rmpath(shadowDir);
if exist(shadowDir, 'dir')
    rmdir(shadowDir, 's');
end
rehash;

function tf = helper_gpuAvailable()
tf = false;
try
    tf = license('test', 'Distrib_Computing_Toolbox') && gpuDeviceCount > 0;
catch
    tf = false;
end
