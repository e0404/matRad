function test_suite = test_generateStf
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
    
    function test_generateStf_photons
        load TG119.mat;
        pln = helper_basicPln('photons',ct,cst);
        
        stf = matRad_generateStf(ct,cst,pln);
        assertTrue(isstruct(stf));
        assertEqual(numel(stf),numel(pln.propStf.gantryAngles));

    function test_generateStf_protons
        load TG119.mat;
        pln = helper_basicPln('protons',ct,cst);
        
        stf = matRad_generateStf(ct,cst,pln);
        assertTrue(isstruct(stf));
        assertEqual(numel(stf),numel(pln.propStf.gantryAngles));

    function test_generateStf_helium
        load TG119.mat;
        pln = helper_basicPln('helium',ct,cst);
        
        stf = matRad_generateStf(ct,cst,pln);
        assertTrue(isstruct(stf));
        assertEqual(numel(stf),numel(pln.propStf.gantryAngles)); 

    function test_generateStf_carbon
        load TG119.mat;
        pln = helper_basicPln('carbon',ct,cst);
        
        stf = matRad_generateStf(ct,cst,pln);
        assertTrue(isstruct(stf));
        assertEqual(numel(stf),numel(pln.propStf.gantryAngles)); 

    function test_generateStf_noTargetObjectives
        load TG119.mat;
        [cst{:,6}] = deal([]);
        pln = helper_basicPln('photons',ct,cst);
                
        stf = matRad_generateStf(ct,cst,pln);
        assertTrue(isstruct(stf));
        assertEqual(numel(stf),numel(pln.propStf.gantryAngles));

    
    function pln = helper_basicPln(modality,ct,cst)
        pln.radiationMode = modality;
        pln.numOfFractions = 30;
        pln.machine = 'Generic';
        pln.multScen = matRad_NominalScenario(ct);
        pln.propStf.gantryAngles = 0;
        pln.propStf.couchAngles = 0;
        pln.propStf.bixelWidth = 5;
        pln.propStf.isoCenter = matRad_getIsoCenter(cst,ct);
        

        