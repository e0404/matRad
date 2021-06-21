classdef matRad_generateBrachyStfTest < matlab.unittest.TestCase
    properties
        OriginalPath
    end
    
    methods (TestMethodSetup)
        function addBrachyToPath(testCase)
            testCase.OriginalPath = path;
            addpath(fileparts(pwd));
            addpath(fullfile(fileparts(fileparts(pwd)),'phantoms'));
        end
    end
    
    methods (TestMethodTeardown)
        function restorePath(testCase)
            path(testCase.OriginalPath);
        end
    end
        
    methods(Test)
        
        
        function rightSeedPositions(testCase)
            % the fucntion should produce the same seed positions as in
            % this simple precomputed example
            load PROSTATE.mat ct cst;
            load plnSimpleTestStruct.mat pln;
            
            stf = matRad_generateBrachyStf(ct, cst, pln);
            actSolutionX = stf.seedPosX;
            actSolutionY = stf.seedPosY;
            actSolutionZ = stf.seedPosZ;
            expSolutionX = [0, 0; 0, 0; 1, 1; 1, 1];
            expSolutionY = [-100, -100; -99, -99; -100, -100; -99, -99];
            expSolutionZ = [1, 2; 1, 2; 1, 2; 1, 2];
            testCase.verifyEqual(actSolutionX,expSolutionX)
            testCase.verifyEqual(actSolutionY,expSolutionY)
            testCase.verifyEqual(actSolutionZ,expSolutionZ)
        end
                
        function outOfBounds(testCase)
            % function should throw an error if it produces seeds out of
            % the ct bounds
            load PROSTATE.mat ct cst;
            load plnSimpleTestStruct.mat pln;
            pln.propStf.orientation.offset = [500,500,500];
            [~] = verifyError(testCase,@()matRad_generateBrachyStf(ct, cst, pln), ...
                'matRad:Error');
        end
         
    end
    
end