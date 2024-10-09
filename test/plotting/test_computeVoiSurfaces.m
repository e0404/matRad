function test_suite = test_computeVoiSurfaces
    %The output should always be test_suite, and the function name the same as
    %your file name
       
    %To collect all tests defined below, this is needed in newer Matlab
    %versions. test_functions will collect function handles to below test
    %functions
    test_functions=localfunctions(); 
    
    % This will initialize the test suite, i.e., take the functions from
    % test_functions, check if they contain "test", convert them into a MOxUnit
    % Test Case, and add them to the test-runner
    initTestSuite;
    
 
    function test_compute

        testpatient = load('photons_testData.mat');
        cst = matRad_computeAllVoiSurfaces(testpatient.ct,testpatient.cst)
        assertTrue(isa(cst,'cell'));
        assertTrue(all(iscell(cst(:,8))))
        assertTrue(all(cellfun(@numel,cst(:,8)) == 2)) % Normals and Vertices        
        assertTrue(isequal(cst(:,1:6),testpatient.cst(:,1:6)))