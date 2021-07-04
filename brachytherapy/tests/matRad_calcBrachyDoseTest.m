classdef matRad_calcBrachyDoseTest < matlab.unittest.TestCase
    % tests all the functions called during dose calculation
     properties
        OriginalPath
    end
    
    methods (TestMethodSetup)
        function addBrachyToPath(testCase)
            testCase.OriginalPath = path;
            addpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))); %matRad
            addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'phantoms')); %matRad/phantoms
        end
    end
    
    methods (TestMethodTeardown)
        function restorePath(testCase)
            path(testCase.OriginalPath);
        end
    end
    
    methods(Test)
        % test if dij struct has required shape
        function rightSructure(testCase)
            load PROSTATE.mat ct cst;
            load calcDoseSimpleTestStruct.mat pln stf;
            pln.propDoseCalc.durationImplanted = Inf;
            
       
            dij = matRad_calcBrachyDose(ct,stf,pln,cst);
            testCase.verifyTrue(isfield(dij,'doseGrid'));
            testCase.verifyTrue(isfield(dij,'physicalDose'));
            testCase.verifyTrue(isfield(dij,'totalNumOfBixels'));
            testCase.verifyTrue(iscell(dij.physicalDose));
            
            clear ct cst pln;
        end
        
        function getDistMtx(testCase)
            % check calculated distance matrix against values calculated
            % on a spreadsheet (same folder as test) for 3x3 test point mtx
            seedPoints.x = [2,1,2];
            seedPoints.y = [2,3,1];
            seedPoints.z = [1,5,1];
            dosePoints.x = [0,1,1];
            dosePoints.y = [0,2,3];
            dosePoints.z = [0,1,2];
            verifyDistanceMatrix = [3,5.9,2.4;1,4.1,1.4;1.7,3,2.45];           
            [funcDistanceMatrix,~] = getDistanceMatrix(seedPoints,dosePoints);
            
            testCase.verifyEqual(funcDistanceMatrix.dist,verifyDistanceMatrix,'AbsTol',0.1)
        end
        
        function getThetMtx(testCase)
            % check calculated angle matrix against values calculated
            % on a spreadsheet (same folder as test) for 3x3 test point mtx
            seedPoints.x = [2,1,2];
            seedPoints.y = [2,3,1];
            seedPoints.z = [1,5,1];
            dosePoints.x = [0,1,1];
            dosePoints.y = [0,2,3];
            dosePoints.z = [0,1,2];
            SeedDirection = sqrt(1/3)*ones(1,3);
            verifyThetaMatrix = [164.2,151.4,160.5;125.3,134.4,90;70.5,125.3,61.9];
            [DistanceMatrix,~] = getDistanceMatrix(seedPoints,dosePoints);
            [funcThetaMatrix,~] = getThetaMatrix(SeedDirection,DistanceMatrix);
            
            testCase.verifyEqual(funcThetaMatrix,verifyThetaMatrix,'AbsTol',0.1)
        end
         
        function plausibleDoseRate1D(testCase)
            % verify if the calculated dose fits precalculated TG43
            % formalism (spreadsheet in test folder)
            load brachy_LDRTest machine;
            
            
        end
        
    end
end
