function test_suite = test_ignoreOutsideDensities

test_functions = localfunctions();

initTestSuite;

function test_higherEnergiesInStf
testData = load('protons_testData.mat');

testData.ct.cubeHU{1}(1,:,:) = -500;

stf_noIgnore = matRad_generateStf(testData.ct,testData.cst,testData.pln);

testData.pln.propStf.ignoreOutsideDensities = true;

stf_Ignore = matRad_generateStf(testData.ct,testData.cst,testData.pln);

% check if energies for beam 1 are higher
for i = 1:size(stf_noIgnore(1).ray,2)
    assertTrue( stf_Ignore(1).ray(i).energy(1)<stf_noIgnore(1).ray(i).energy(1));
end

for i = 1:size(stf_noIgnore(2).ray,2)
    assertTrue( isequal(stf_Ignore(2).ray(i).energy,stf_noIgnore(2).ray(i).energy));
end
