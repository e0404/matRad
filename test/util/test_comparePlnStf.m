function test_suite = test_comparePlnStf

test_functions=localfunctions();

initTestSuite;

function [pln,stf] = helper_getDummyBasicPlnStf()
    pln = struct('radiationMode', 'photons',...
        'propStf', struct('numOfBeams', 2, 'gantryAngles', [0, 90],'couchAngles', [0, 90],'bixelWidth', 5, 'isoCenter', [0, 0, 0; 1, 1, 1]));
    stf = [struct('gantryAngle', 0,'couchAngle', 0,'bixelWidth', 5, 'isoCenter', [0, 0, 0], 'radiationMode', 'photons'),...
        struct('gantryAngle', 90,'couchAngle', 90,'bixelWidth', 5, 'isoCenter', [1, 1, 1], 'radiationMode', 'photons')];

        % Match should not work if pln does not have propStf
function test_missingPropStf
    pln = struct();
    stf = struct();
    [allMatch, msg] = matRad_comparePlnStf(pln, stf);
    assertFalse(allMatch);
    assertFalse(isempty(msg));
    assertTrue(ischar(msg));

% Test case for matching gantry angles
function test_matching
    [pln, stf] = helper_getDummyBasicPlnStf();    
    [allMatch, msg] = matRad_comparePlnStf(pln, stf);
    assertTrue(allMatch);
    assertTrue(isempty(msg));

% Test case for non-matching gantry angles
function test_nonMatchingGantryAngles
    [pln, stf] = helper_getDummyBasicPlnStf();
    pln.propStf.gantryAngles = [0, 75];
    [allMatch, msg] = matRad_comparePlnStf(pln, stf);
    assertFalse(allMatch);
    assertFalse(isempty(msg));
    assertTrue(ischar(msg));

% Test case for non-matching couch angles
function test_nonMatchingCouchAngles
    [pln, stf] = helper_getDummyBasicPlnStf();
    pln.propStf.couchAngles = [0, 75];
    [allMatch, msg] = matRad_comparePlnStf(pln, stf);
    assertFalse(allMatch);
    assertFalse(isempty(msg));
    assertTrue(ischar(msg));

%Test case for wrong number of beams provided
function test_wrongNumberOfAngles
    [pln, stf] = helper_getDummyBasicPlnStf();
    pln.propStf.numOfBeams = 3;
    [allMatch, msg] = matRad_comparePlnStf(pln, stf);
    assertFalse(allMatch);
    assertFalse(isempty(msg));
    assertTrue(ischar(msg));

% Test case for non-matching bixel width
function test_nonMatchingBixelWidth
    [pln, stf] = helper_getDummyBasicPlnStf();
    pln.propStf.bixelWidth = 10;
    [allMatch, msg] = matRad_comparePlnStf(pln, stf);
    assertFalse(allMatch);
    assertFalse(isempty(msg));
    assertTrue(ischar(msg));

% Test case for non-matching radiation mode
function test_nonMatchingRadiationMode
    [pln, stf] = helper_getDummyBasicPlnStf();
    pln.radiationMode = 'protons';
    [allMatch, msg] = matRad_comparePlnStf(pln, stf);
    assertFalse(allMatch);
    assertFalse(isempty(msg));
    assertTrue(ischar(msg));

% Test case for non-matching isocenters
function test_nonMatchingIsocenters
    [pln, stf] = helper_getDummyBasicPlnStf();
    pln.propStf.isoCenter = [0, 0, 0; 1, 1, 2];
    [allMatch, msg] = matRad_comparePlnStf(pln, stf);
    assertFalse(allMatch);
    assertFalse(isempty(msg));
    assertTrue(ischar(msg));
