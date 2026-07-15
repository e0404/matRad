function test_suite = test_ignoreOutsideDensities

test_functions = localfunctions();

initTestSuite;

function test_lowerEnergiesInStf
testData = load('protons_testData.mat');

testData.ct.cubeHU{1}(1,:,:) = -500;

stf_noIgnore = matRad_generateStf(testData.ct,testData.cst,testData.pln);

testData.pln.propStf.ignoreOutsideDensities = true;

stf_Ignore = matRad_generateStf(testData.ct,testData.cst,testData.pln);

% ignoring the outside density slab reduces the radiological depth for
% beam 1 (entering through the modified face), so its energies must be lower
for i = 1:size(stf_noIgnore(1).ray,2)
    assertTrue( stf_Ignore(1).ray(i).energy(1)<stf_noIgnore(1).ray(i).energy(1));
end

% for beam 2 the slab is distal to the target, so the energies are unaffected
for i = 1:size(stf_noIgnore(2).ray,2)
    assertTrue( isequal(stf_Ignore(2).ray(i).energy,stf_noIgnore(2).ray(i).energy));
end
