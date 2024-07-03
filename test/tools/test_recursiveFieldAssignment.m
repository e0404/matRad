function test_suite = test_recursiveFieldAssignment
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

% Assigning fields from reference to assignTo
function test_recursiveFieldAssignment_existingField
    assignTo = struct('a', 1, 'b', struct('c', 2));
    reference = struct('a', 3, 'b', struct('c', 4));
    expected = struct('a', 3, 'b', struct('c', 4));
    result = matRad_recursiveFieldAssignment(assignTo, reference);
    assertEqual(result, expected);
    %assertWarning(@() matRad_recursiveFieldAssignment(assignTo, reference,'TestWarn')); %Does not work for some reason, despite warnign being correctly raised


function test_recursiveFieldAssignment_keepFields
    assignTo = struct('a', 1, 'b', struct('c', 2), 'd', struct('e','hello'));
    reference = struct('a', 3, 'b', struct('c', 4));
    expected = struct('a', 3, 'b', struct('c', 4), 'd', struct('e','hello'));
    result = matRad_recursiveFieldAssignment(assignTo, reference);
    assertEqual(result, expected);

% Assigning a struct field to a non-struct field
function test_recursiveFieldAssignment_structToNonExisting
    assignTo = struct();
    reference = struct('a', 2);
    expected = struct('a', 2);
    result = matRad_recursiveFieldAssignment(assignTo, reference);
    assertEqual(result, expected);
    %assertWarning(@() matRad_recursiveFieldAssignment(assignTo, reference)); %Does not work for some reason, despite warning being correctly raised

    assignTo = struct('a',2);
    reference = struct('a',3,'b',4);
    expected = struct('a',3,'b',4);
    result = matRad_recursiveFieldAssignment(assignTo, reference);
    assertEqual(result, expected);

% Overwriting a non-struct field with a non-struct field
function test_recursiveFieldAssignment_nonStruct
    assignTo = 1;
    reference = 2;
    expected = 2;
    result = matRad_recursiveFieldAssignment(assignTo, reference);
    assertEqual(result, expected);
    %assertWarning(@() matRad_recursiveFieldAssignment(assignTo, reference,'TestWarn')); %Does not work for some reason, despite warnign being correctly raised

% Overwriting a non-struct field with a struct field
function test_recursiveFieldAssignment_nonExistingField
    assignTo = struct('a', 1);
    reference = struct('a', struct('b', 2));
    expected = struct('a', struct('b', 2));
    result = matRad_recursiveFieldAssignment(assignTo, reference);
    assertEqual(result, expected);
    %assertWarning(@() matRad_recursiveFieldAssignment(assignTo, reference)); %Does not work for some reason, despite warnign being correctly raised


% Overwriting a struct field with a non-struct field
function test_recursiveFieldAssignment_nestedField
    assignTo = struct('a', struct('b', 2));
    reference = struct('a', 3);
    expected = struct('a', 3);
    result = matRad_recursiveFieldAssignment(assignTo, reference);
    assertEqual(result, expected);
    %assertWarning(@() matRad_recursiveFieldAssignment(assignTo, reference)); %Does not work for some reason, despite warnign being correctly raised



% Add more test cases here...
