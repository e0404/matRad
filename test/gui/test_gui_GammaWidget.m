function test_suite = test_gui_GammaWidget
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

function test_GammaWidget_constructor
    % Test the constructor of matRad_GammaWidget class
    h = matRad_GammaWidget();
    assertTrue(isa(h, 'matRad_GammaWidget'));
    assertTrue(isa(h, 'matRad_Widget'));
    delete(h);

    p = uipanel();
    h = matRad_GammaWidget(p);
    assertTrue(isa(h, 'matRad_GammaWidget'));
    assertTrue(isa(h, 'matRad_Widget'));
    delete(h);
    close(get(p,'Parent'));

%TODO: Test Buttons / visibility depending on data