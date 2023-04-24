classdef matRad_generateBrachyStfTest < matlab.unittest.TestCase
    properties
        OriginalPath
        plnTest
    end
    
    methods (TestMethodSetup)
        function addBrachyToPath(testCase)
            testCase.OriginalPath = path;
            addpath(fileparts(pwd));
            addpath(fullfile(fileparts(fileparts(pwd)),'phantoms'));
            
            % set up an example pln with relevant fields 
            pln.radiationMode   = 'brachy';
            pln.machine         = 'HDR';
            % geometry settings
            load PROSTATE.mat ct cst;
            pln.propStf.templateRoot             = matRad_getTemplateRoot(ct,cst);
            pln.propStf.needle.seedDistance      = 1; % [mm] seed distance on needle
            pln.propStf.needle.seedsNo           = 2; % number of seeds per needle
            pln.propStf.template.numOfXPoints  = 2;
            pln.propStf.template.numOfYPoints = 2;
            pln.propStf.template.xScale        = 1; % [mm] distance of neighbouring points
            pln.propStf.template.yScale       = 1; % [mm] distance of neighbouring points
            pln.propStf.template.activeNeedles = [0 0 0 0 0 0 0 0 0 0 0 0 0;... % 7.0
                                                  0 0 0 0 0 0 0 0 0 0 0 0 0;... % 6.5
                                                  0 0 0 0 1 1 0 1 1 0 0 0 0;... % 6.0
                                                  1 0 1 0 1 0 0 0 1 0 1 0 1;... % 5.5
                                                  0 1 0 1 0 1 0 1 0 1 0 1 0;... % 5.0
                                                  1 0 1 0 1 0 0 0 1 0 1 0 1;... % 4.5
                                                  0 1 0 1 0 1 0 1 0 1 0 1 0;... % 4.0
                                                  1 0 1 0 1 0 1 0 1 0 1 0 1;... % 4.5
                                                  0 1 0 1 0 1 0 1 0 1 0 1 0;... % 3.0
                                                  0 0 1 0 1 0 0 0 1 0 1 0 0;... % 2.5
                                                  0 0 0 1 0 0 0 0 0 1 0 0 0;... % 2.0
                                                  0 0 0 0 0 0 0 0 0 0 0 0 0;... % 1.5                                                                          
                                                  0 0 0 0 0 0 0 0 0 0 0 0 0];   % 1.0
                                                 %A a B b C c D d E e F f G
            pln.propStf.bixelWidth      = 10;
            testCase.plnTest = pln;         
        end
    end
    
    methods (TestMethodTeardown)
        function restorePath(testCase)
            path(testCase.OriginalPath);
        end
    end
        
    methods(Test)
        
        function rightOutput(testCase)
            load PROSTATE.mat ct cst;
            stf = matRad_generateBrachyStf(ct,cst,testCase.plnTest,0);
            testCase.verifyTrue(isfield(stf,'radiationMode'));
            testCase.verifyTrue(isfield(stf,'numOfSeedsPerNeedle'));
            testCase.verifyTrue(isfield(stf,'numOfNeedles'));
            testCase.verifyTrue(isfield(stf,'totalNumOfBixels'));
            testCase.verifyTrue(isfield(stf,'template'));
            testCase.verifyTrue(isfield(stf,'seedPoints'));
        end

            
        
         
    end
    
end