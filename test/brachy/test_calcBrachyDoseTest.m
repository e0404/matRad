function test_suite = test_calcBrachyDoseTest

test_functions = localfunctions();

initTestSuite;



% test if dij struct has required shape
function rightOutput(testCase)
load PROSTATE.mat ct cst;
load examplePlnAndStf.mat pln stf;
pln.propDoseCalc.durationImplanted = Inf;
            
       
dij = testCase.engine.setCalcDose(ct,cst,stf);
testCase.verifyTrue(isfield(dij,'doseGrid'));
testCase.verifyTrue(isfield(dij,'physicalDose'));
% testCase.verifyTrue(isfield(dij,'basedata'));
testCase.verifyTrue(iscell(dij.physicalDose));
            
clear ct cst pln;
end

