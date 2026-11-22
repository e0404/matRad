function test_suite = test_isgpuarrayCompat

test_functions = localfunctions();

initTestSuite;

function test_isgpuarrayCompatHostArrays
assertFalse(matRad_isgpuarrayCompat(1:3));
assertFalse(matRad_isgpuarrayCompat(single(ones(2))));
assertFalse(matRad_isgpuarrayCompat(true));
assertFalse(matRad_isgpuarrayCompat({1, 2}));
assertFalse(matRad_isgpuarrayCompat(struct('a', 1)));

function test_isgpuarrayCompatGpuArray
matRad_cfg = MatRad_Config.instance();
if ~matRad_cfg.isMatlab || exist('gpuDeviceCount', 'file') == 0 || gpuDeviceCount() == 0
    moxunit_throw_test_skipped_exception('No GPU available');
end
assertTrue(matRad_isgpuarrayCompat(gpuArray(1:3)));

function test_interp1WithoutIsgpuarray
shadowDir = helper_temporaryFolder('matRad_test_isgpuarrayShadow');
fid = fopen(fullfile(shadowDir, 'isgpuarray.m'), 'w');
fprintf(fid, 'function tf = isgpuarray(~)\nerror(''isgpuarray must not be called'');\nend\n');
fclose(fid);
addpath(shadowDir);
cleanup = onCleanup(@() rmpath(shadowDir));

y = matRad_interp1([0; 1; 2], [0; 10; 20], [0.5; 1.5]);
assertElementsAlmostEqual(y, [5; 15]);
