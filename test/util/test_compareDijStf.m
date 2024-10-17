function test_suite = test_matRad_compareDijStf
test_functions=localfunctions();

initTestSuite;

function [dij,stf] = helper_getDummyBasicDijStf()
    dij = struct('numOfRaysPerBeam', [100, 200], 'numOfBeams', 2);
    stf = [struct('numOfRays',100,'gantryAngle',0), struct('numOfRays',200,'gantryAngle',90)];

% Test case 1: Matching dij and stf
function test_matchingDijStf
    [dij, stf] = helper_getDummyBasicDijStf();
    [allMatch, msg] = matRad_compareDijStf(dij, stf);
    assertTrue(allMatch);
    assertTrue(isempty(msg));
    
% Test case 2: Different number of rays per beam
function test_differentRaysPerBeam
    [dij, stf] = helper_getDummyBasicDijStf();
    dij.numOfRaysPerBeam = [100, 300];
    [allMatch, msg] = matRad_compareDijStf(dij, stf);
    assertFalse(allMatch);
    assertFalse(isempty(msg));
    assertTrue(ischar(msg));

% Test case 3: Different number of beams
function test_differentNumOfBeams
    [dij, stf] = helper_getDummyBasicDijStf();
    dij.numOfBeams = 3;
    [allMatch, msg] = matRad_compareDijStf(dij, stf);
    assertFalse(allMatch);
    assertFalse(isempty(msg));
    assertTrue(ischar(msg));

