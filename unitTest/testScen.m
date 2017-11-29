classdef testScen < matlab.unittest.TestCase
    properties
        % standard allowed relative deviation
        rEps = 1e-5;
        aEpsRough = 1e-2;
        
        % testscenarios
        testCaseFile = 'testCases.mat';
        testCaseHash = 'f3f57dad5fc4676e29d4b9a9b5ff31f3';
    end
    
    methods (Test)
        %% test hash
        function validateDataset(testScen)
            load(testScen.testCaseFile);
            actualHash = DataHash(testCase);
            testScen.fatalAssertEqual(actualHash, testScen.testCaseHash);
        end
        %% test multScen
        function testProbCalcHard(testScen)
            expectedProb = 1;
            actualProb = matRad_calcScenProb([0 0 0 0 0], [0 0 0 0 0], zeros(1,5),'probBins','normDist');
            testScen.assertEqual(actualProb, expectedProb);
            actualProb = matRad_calcScenProb([0 0 0 0 0], [1 1 1 1 1], zeros(1,5),'pointwise','normDist');
            testScen.assertEqual(actualProb, expectedProb);
            
            % simple relative range shift scenario
            rangeShiftScenario = zeros(3,5);
            rangeShiftScenario(1,5) = -0.03;
            rangeShiftScenario(3,5) = 0.03;
            
            expectedProb = [0 1 0]';
            actualProb = matRad_calcScenProb([0 0 0 0 0], [0 0 0 0 0], rangeShiftScenario,'probBins','normDist');
            testScen.assertEqual(actualProb, expectedProb);
            
            expectedProb = [0.2790 0.4420 0.2790]';
            actualProb = matRad_calcScenProb([0 0 0 0 0], [0 0 0 0 0.03], rangeShiftScenario,'probBins','normDist');
            testScen.assertEqual(actualProb, expectedProb, 'AbsTol', testScen.aEpsRough);
        end
        
        function propagateMultScen(testScen)
            ct.numOfCtScen = 1;
            pln.sampling = true;
            
            load(testScen.testCaseFile);
            testSize = numel(testCase.multScenCases);
            for i = 1:testSize
                multScen = testCase.multScenCases(i).multScen;
                actualPln = matRad_setPlanUncertainties(ct, pln, multScen);
                
                actualProb = actualPln.multScen.scenProb;
                expectedProb = testCase.multScenCases(i).prob;
                testScen.assertEqual(actualProb, expectedProb, 'RelTol', testScen.rEps);
                
                actualScenParam = actualPln.multScen.scenForProb;
                expectedScenParam = testCase.multScenCases(i).scenParam;
                testScen.assertEqual(actualScenParam, expectedScenParam, 'RelTol', testScen.rEps);
            end
        end
        
        
        %% statistical tests
        function percentileBig(testScen)
            testSize = 100;
            sampleSize = 10000;
            numPercentiles = 10;
            
            % big sample size
            for i = 1:testSize
                x = rand([1,sampleSize]);
                percentiles = rand([1,numPercentiles]);
                expected = quantile(x,percentiles);
                actual = matRad_weightedQuantile(x,percentiles,[],false,'nearest');
                testScen.assertEqual(actual, expected, 'RelTol', testScen.rEps);
            end
        end
        
        function percentileSmall(testScen)
            testSize = 100;
            sampleSize = 5;
            numPercentiles = 10;

            % small sample size
            for i = 1:testSize
                x = rand([1,sampleSize]);
                percentiles = rand([1,numPercentiles]);
                expected = quantile(x,percentiles);
                actual = matRad_weightedQuantile(x,percentiles,[],false,'nearest');
                testScen.assertEqual(actual, expected, 'RelTol', testScen.rEps);
            end
        end
    end
    
end