classdef matRadBrachy_generateStfTest < matlab.unittest.TestCase
    methods(Test)
        function trivialExample(testCase)
            load PROSTATE.mat ct cst;
            pln.radiationMode            = 'brachy';
            pln.machine                  = 'Generic';
            pln.propStf.template.numOfHorPoints  = 2;
            pln.propStf.template.numOfVertPoints = 2;
            pln.propStf.template.Xscale          = 1; % [mm]
            pln.propStf.template.Yscale          = 1; % [mm]
            pln.propStf.needle.seedDistance      = 1; % [mm]
            pln.propStf.needle.seedsNo           = 2; 
            %unit vectors of displaced, rotated template coordinate system
            pln.propStf.orientation.Xdir    = normalize([1,0,0],'norm');
            pln.propStf.orientation.Ydir    = normalize([0,1,0],'norm');
            pln.propStf.orientation.Zdir    = cross(pln.propStf.orientation.Xdir,pln.propStf.orientation.Ydir);
            pln.propStf.orientation.offset  = [0,0,0]; % [mm]
            
            stf = matRadBrachy_generateStf(ct, cst, pln);
            actSolutionX = stf.seedPosX;
            actSolutionY = stf.seedPosY;
            actSolutionZ = stf.seedPosZ;
            expSolutionX = [0, 0; 0, 0; 1, 1; 1, 1];
            expSolutionY = [0, 0; 1, 1; 0, 0; 1, 1];
            expSolutionZ = [1, 2; 1, 2; 1, 2; 1, 2];
            testCase.verifyEqual(actSolutionX,expSolutionX)
            testCase.verifyEqual(actSolutionY,expSolutionY)
            testCase.verifyEqual(actSolutionZ,expSolutionZ)
        end
        function advancedExample(testCase)
            load PROSTATE.mat ct cst;
            pln.radiationMode            = 'brachy';
            pln.machine                  = 'Generic';
            pln.propStf.template.numOfHorPoints  = 6;
            pln.propStf.template.numOfVertPoints = 8;
            pln.propStf.template.Xscale          = 5; % [mm]
            pln.propStf.template.Yscale          = 5; % [mm]
            pln.propStf.needle.seedDistance      = 20; % [mm]
            pln.propStf.needle.seedsNo           = 10; 
            %unit vectors of displaced, rotated template coordinate system
            pln.propStf.orientation.Xdir = normalize([10,0,1],'norm');
            pln.propStf.orientation.Ydir = normalize([0.02,1,-0.2],'norm');
            pln.propStf.orientation.Zdir = cross(pln.propStf.orientation.Xdir,pln.propStf.orientation.Ydir);
            pln.propStf.orientation.offset = [0,0,0]; % [mm]
            
            stf = matRadBrachy_generateStf(ct, cst, pln);
            actSolutionX = stf.seedPosX;
            actSolutionY = stf.seedPosY;
            actSolutionZ = stf.seedPosZ;
            load validSTF.mat stfEXP
            expSolutionX = stfEXP.seedPosX;
            expSolutionY = stfEXP.seedPosY;
            expSolutionZ = stfEXP.seedPosZ;
            testCase.verifyEqual(actSolutionX,expSolutionX)
            testCase.verifyEqual(actSolutionY,expSolutionY)
            testCase.verifyEqual(actSolutionZ,expSolutionZ)
        end
        
        
        
%         function nonnumericInput(testCase)
%             testCase.verifyError(@()quadraticSolver(1,'-3',2), ...
%                 'quadraticSolver:InputMustBeNumeric')
%         end
    end
end