function test_suite = test_io_nrrd
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

function test_nrrd_basic_writeReadCube
    load TG119.mat
    tmpDir = tempdir();
    tmpDir = [tmpDir filesep 'matRad_test_nrrd' filesep];
    if ~exist(tmpDir,'dir')
        mkdir(tmpDir);
    end

    tmpFilePath = [tmpDir filesep 'TG119.nrrd'];
    
    datatype = 'double';

    metadata = struct;
    metadata.resolution = ct.resolution;
    
    %Check the write
    matRad_writeCube(tmpFilePath,ct.cubeHU{1},datatype,metadata);
    assertEqual(exist(tmpFilePath,'file'),2);
    
    %check the read
    [readCube,readMeta] = matRad_readCube(tmpFilePath);

    assertElementsAlmostEqual(ct.cubeHU{1},readCube);    
    assertEqual(readMeta.resolution,[ct.resolution.x ct.resolution.y ct.resolution.z]);
    assertEqual(readMeta.cubeDim,ct.cubeDim);
    

