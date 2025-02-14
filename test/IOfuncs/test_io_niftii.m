function test_suite = test_io_niftii
%The output should always be test_suite, and the function name the same as
%your file name

%% Header
% The header is required to be in this format for automatic test collection
% by MOxUnit

%To collect all tests defined below, this is needed in newer Matlab
%versions. test_functions will collect function handles to below test
%functions
test_functions=localfunctions(); 

% This will initialize the test suite, i.e., take the functions from
% test_functions, check if they contain "test", convert them into a MOxUnit
% Test Case, and add them to the test-runner
initTestSuite;

%% Custom Tests
% Tests use assert*-like Functions to check outputs etc:
% assertTrue(a) - a is true
% assertFalse(a) - a is false
% assertEqual(a,b) - a and be are equal (isequal)
% assertElementsAlmostEqual(a,b) - numerical test for all vector / matrix elements. Has Additional arguments for absolute / relative tolerance 
% assertVectorsAlmostEqual(a,b) - numerical test using vector norm
% assertExceptionThrown(f,id) - test if exception of id is thrown. Take care of Octave issues with exception id (or don't provide id)
% Check MOxUnit for more information or look at other tests

function test_niftii_basic_writeReadCube
    if moxunit_util_platform_is_octave
        moxunit_throw_test_skipped_exception('niftiread not available for Octave!');
    end

    load TG119.mat

    tmpDir = helper_temporaryFolder('matRad_test_nii');
    tmpFilePath = fullfile(tmpDir,'TG119.nii');
    
    datatype = 'double';

    metadata = struct;
    metadata.resolution = ct.resolution;
    metadata.compress = false;
    
    %Check the write
    matRad_writeCube(tmpFilePath,ct.cubeHU{1},datatype,metadata);
    assertEqual(exist(tmpFilePath,'file'),2);
    
    %check the read
    [readCube,readMeta] = matRad_readCube(tmpFilePath);

    assertElementsAlmostEqual(ct.cubeHU{1},readCube);    
    assertEqual(readMeta.resolution,[ct.resolution.x ct.resolution.y ct.resolution.z]);
    assertEqual(readMeta.cubeDim,ct.cubeDim);

function test_niftii_compressed_writeReadCube
    if moxunit_util_platform_is_octave
        moxunit_throw_test_skipped_exception('niftiread not available for Octave!');
    end

    load TG119.mat
    
    tmpDir = helper_temporaryFolder('matRad_test_nii');
    tmpFilePath = fullfile(tmpDir,'TG119_compressed.nii');
    
    datatype = 'double';

    metadata = struct;
    metadata.resolution = ct.resolution;
    metadata.compress = true;
    
    %Check the write
    matRad_writeCube(tmpFilePath,ct.cubeHU{1},datatype,metadata);

    assertEqual(exist([tmpFilePath '.gz'],'file'),2);
    
    %check the read
    [readCube,readMeta] = matRad_readCube([tmpFilePath '.gz']);

    assertElementsAlmostEqual(ct.cubeHU{1},readCube);    
    assertEqual(readMeta.resolution,[ct.resolution.x ct.resolution.y ct.resolution.z]);
    assertEqual(readMeta.cubeDim,ct.cubeDim);
    

