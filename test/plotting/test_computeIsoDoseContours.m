function test_suite = test_computeIsoDoseContours
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
        dim = [10,13,15];
        dose_cube = rand(dim);
        levels = [0.1,0.5,0.95];
        isoDoseContours = matRad_computeIsoDoseContours(dose_cube,levels);
        assertTrue(iscell(isoDoseContours));
        assertEqual(size(isoDoseContours),[max(dim),3])
        sz = cellfun(@(c) size(c,1),isoDoseContours(:));
        assertTrue(all(sz == 0 | sz == 2))
        