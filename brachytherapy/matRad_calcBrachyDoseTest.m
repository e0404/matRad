classdef matRad_calcBrachyDoseTest < matlab.unittest.TestCase
    methods(Test)
        function rightSructure(testCase)
            load PROSTATE.mat ct cst;
            % the following is an arbitraty minimal example for a brachy
            % pln struct
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
            
            dij = matRad_calcBrachyDose(ct,stf,pln,cst);
            testCase.verifyTrue(isfield(dij,'doseGrid'));
            testCase.verifyTrue(isfield(dij,'physicalDose'));
            testCase.verifyTrue(isfield(dij,'totalNumOfBixels'));
            testCase.verifyTrue(iscell(dij.physicalDose));
        end
        

    end
end
