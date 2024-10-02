function test_suite = test_DicomImporter
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

function test_DicomImporter_classattachment
    path = 'C:\Users\r114m\Documents\GitHub\matRad\userdata\dicomExport\HEAD_AND_NECK';
    h = matRad_DicomImporter(path);
    assertTrue(isa(h,'matRad_DicomImporter'));
    assertTrue(isa(h,'handle'));

function test_DicomImporter_loadFiles
    path = 'C:\Users\r114m\Documents\GitHub\matRad\userdata\dicomExport\HEAD_AND_NECK';
    h = matRad_DicomImporter(path);
    assertTrue(isequal(h.patDir, path));
    assertTrue(~isempty(h.allfiles));

    NumberOfCtFiles = numel(nonzeros(strcmp(h.allfiles(:,2),'CT')));
    NumberOfRtssFiles = numel(nonzeros(strcmpi(h.allfiles(:,2),'rtstruct'))); 
    NumberOfRtPlanFiles = numel(nonzeros(strcmpi(h.allfiles(:,2),'rtplan')));
    NumberOfRtDoseFiles = numel(nonzeros(strcmpi(h.allfiles(:,2),'rtdose')));

    assertTrue(isequal(NumberOfCtFiles, size(h.importFiles.ct, 1)));
    assertTrue(isequal(NumberOfRtssFiles, size(h.importFiles.rtss, 1)));
    assertTrue(isequal(NumberOfRtPlanFiles, size(h.importFiles.rtplan,1)));
    assertTrue(isequal(NumberOfRtDoseFiles, size(h.importFiles.rtdose,1)));
    
    resBool = true;
    if isempty(h.importFiles.resx) || isempty(h.importFiles.resy) || isempty(h.importFiles.resz)
       resBool = false;
    end
    assertTrue(resBool);



    
    





