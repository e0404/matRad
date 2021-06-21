% This is a script to automatically run all tests the test - folder.
%% setup
if ~exist('matRad_cfg','var')
addpath(fileparts(pwd));
matRad_rc;
end

%% run test for stf
testStf = matRad_generateBrachyStfTest;
testResultStf = run(testStf)

%% run test for CalsDose
testCalcDose = matRad_calcBrachyDoseTest;
testResultCalcDose = run(testCalcDose)