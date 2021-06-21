classdef matRad_getDoseRate2DTest < matlab.unittest.TestCase
%     methods(TestMethodSetup)
%         function createFigure(testCase)
%             comment
%             testCase.TestFigure = figure;
%         end
%     end
%  
%     methods(TestMethodTeardown)
%         function closeFigure(testCase)
%             close(testCase.TestFigure)
%         end
%     end
    
    methods(Test)
        function ValuePlausibility1(testCase)
            % Compares the function result 'actDose' to a 
            % manually calculated one 'calcDose'.
            load brachy_Generic machine;
            sourceObj = Source(machine.data);
            actDose = sourceObj.getDoseRate2D(3,60);
            calcDose = 0.1310;
            tolerance = calcDose*0.1;
            testCase.verifyTrue((actDose-calcDose)<tolerance);
        end
        function ValuePlausibility2(testCase)
            % Compares the function result 'actDose' to a 
            % manually calculated one 'calcDose'.
            load brachy_LDR machine;
            sourceObj = Source(machine.data);
            actDose = sourceObj.getDoseRate2D(5,deg2rad(10));
            calcDose = 0.0086;
            tolerance = calcDose*0.1;
            testCase.verifyTrue((actDose-calcDose)<tolerance);
        end
    end
end
