classdef matRad_generateStfTest < matlab.unittest.TestCase
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
            
            stf = matRad_generateStf(ct, cst, pln);
            actSolutionX = stf.seedPoints.x;
            actSolutionY = stf.seedPoints.y;
            actSolutionZ = stf.seedPoints.z;
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
            pln.propStf.shiftRotMtx(1:3,3) = [500,500,500];
            [~] = verifyError(testCase,@()matRad_generateStf(ct, cst, pln), ...
                'matRad:Error');
        end
         
    end
    
end