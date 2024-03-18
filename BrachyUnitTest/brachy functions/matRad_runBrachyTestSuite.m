% This is a script to automatically run all tests this test - folder.
%% setup
if ~exist('matRad_cfg','var')
addpath(fileparts(fileparts(pwd)));
matRad_cfg =  MatRad_Config.instance();
matRad_rc;
end

%% run test for stf
% this test checks if all relevant fields are contained in the output.

testStf = matRad_generateBrachyStfTest;
testResultStf = run(testStf)

%% run test for CalcDose
% this test suite does five seperate tests:
% 1. Test if the output of matRad_calcBrachyDose contains all fields
% 2. Test if matRad_getDistanceMatrix results conform to hand calculated
% results
% 3. Test if matRad_getThetaMatrix results conform to hand calculated
% results
% 4. Test if  matRad_getDoseRate1D_poly results conform to precalculated
% results
% 5. Test if  matRad_getDoseRate2D_poly results conform to precalculated
% results

testCalcDose = matRad_calcBrachyDoseTest;
testResultCalcDose = run(testCalcDose)