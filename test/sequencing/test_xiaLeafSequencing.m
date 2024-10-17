function test_suite = test_xiaLeafSequencing
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

function [resultGUI,stf,dij] = helper_getTestData()
    p = load('photons_testData.mat');
    resultGUI = p.resultGUI;
    stf = p.stf;
    dij = p.dij;


function test_run_sequencing_basic
    [resultGUI,stf,dij] = helper_getTestData();
    fn_old = fieldnames(resultGUI);
    
    numOfLevels = [1,10];

    for levels = numOfLevels
        resultGUI_sequenced = matRad_xiaLeafSequencing(resultGUI,stf,dij,levels);
        
        fn_new = fieldnames(resultGUI_sequenced);
        for i = 1:numel(fn_old)
            assertTrue(any(strcmp(fn_old{i},fn_new)));
            assertEqual(resultGUI.(fn_old{i}),resultGUI_sequenced.(fn_old{i}));
        end
        
        % Basic additions to resultGUI
        assertTrue(isvector(resultGUI_sequenced.wSequenced));
        assertTrue(isstruct(resultGUI_sequenced.apertureInfo));
        assertTrue(isstruct(resultGUI_sequenced.sequencing));

        %Sequencing Struct
        seq = resultGUI_sequenced.sequencing;
        assertTrue(isstruct(seq.beam));
        assertTrue(numel(seq.beam) == numel(stf));
        for i = 1:numel(seq.beam)
            assertTrue(isscalar(seq.beam(i).numOfShapes));
            assertTrue(isnumeric(seq.beam(i).shapes));
            shapeSize = size(seq.beam(i).shapes);
            assertEqual(shapeSize(3),seq.beam(i).numOfShapes);
            assertTrue(isvector(seq.beam(i).shapesWeight));
            assertTrue(isvector(seq.beam(i).bixelIx));
            assertTrue(ismatrix(seq.beam(i).fluence));
            assertEqual(size(seq.beam(i).fluence),shapeSize([1 2]));
        end
        
        %ApertureInfo Sturct
        apInfo = resultGUI_sequenced.apertureInfo;
        assertTrue(isscalar(apInfo.bixelWidth));
        assertTrue(isscalar(apInfo.numOfMLCLeafPairs));
        assertTrue(isscalar(apInfo.totalNumOfBixels));
        assertTrue(isscalar(apInfo.totalNumOfShapes));
        assertTrue(isscalar(apInfo.totalNumOfLeafPairs));
        assertTrue(isvector(apInfo.apertureVector));
        assertTrue(ismatrix(apInfo.mappingMx));
        assertTrue(ismatrix(apInfo.limMx));
        assertTrue(isstruct(apInfo.beam))
        assertTrue(numel(apInfo.beam) == numel(stf));
    end


    

    





        