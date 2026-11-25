function test_suite = test_resolveClassSelection
% Tests for matRad_resolveClassSelection, which selects a class from a list
% by its shortName and explains why a requested class cannot be used (#943)

test_functions = localfunctions();

initTestSuite;

function pln = helper_photonPln()
pln = struct('radiationMode', 'photons', 'machine', 'Generic');

function [availableClasses, allClasses] = helper_engineLists(pln)
availableClasses = DoseEngines.matRad_DoseEngineBase.getAvailableEngines(pln);
allClasses = DoseEngines.matRad_DoseEngineBase.getAvailableEngines();

function test_resolveNameIsCaseInsensitive
pln = helper_photonPln();
[availableClasses, allClasses] = helper_engineLists(pln);
[ix, reason] = matRad_resolveClassSelection('svdpb', availableClasses, allClasses, pln);
assertTrue(isscalar(ix));
assertEqual(availableClasses(ix).shortName, 'SVDPB');
assertTrue(isempty(reason));

function test_resolveStringName
if moxunit_util_platform_is_octave()
    moxunit_throw_test_skipped_exception('Octave has no string class');
end
pln = helper_photonPln();
[availableClasses, allClasses] = helper_engineLists(pln);
ix = matRad_resolveClassSelection(string('SVDPB'), availableClasses, allClasses, pln);
assertEqual(availableClasses(ix).shortName, 'SVDPB');

function test_resolveListTakesFirstAvailableName
% a cell array of names, as in the defaults of MatRad_Config, selects the
% first name that is available
pln = helper_photonPln();
[availableClasses, allClasses] = helper_engineLists(pln);
[ix, reason] = matRad_resolveClassSelection({'MCsquare', 'SVDPB', 'HongPB'}, availableClasses, allClasses, pln);
assertEqual(availableClasses(ix).shortName, 'SVDPB');
assertTrue(isempty(reason));

function test_resolveDoesNotEvaluateAllClassesOnMatch
pln = helper_photonPln();
availableClasses = helper_engineLists(pln);
ix = matRad_resolveClassSelection('SVDPB', availableClasses, @() error('test:called', 'must not be called'), pln);
assertEqual(availableClasses(ix).shortName, 'SVDPB');

function test_resolveAcceptsFunctionHandleForAllClasses
pln = helper_photonPln();
availableClasses = helper_engineLists(pln);
[ix, reason] = matRad_resolveClassSelection('MCsquare', availableClasses, @() DoseEngines.matRad_DoseEngineBase.getAvailableEngines(), pln);
assertTrue(isempty(ix));
assertFalse(isempty(strfind(reason, 'photons')));

function test_resolveUnknownNameListsAvailableNames
pln = helper_photonPln();
[availableClasses, allClasses] = helper_engineLists(pln);
[ix, reason] = matRad_resolveClassSelection('MCsqare', availableClasses, allClasses, pln);
assertTrue(isempty(ix));
assertFalse(isempty(strfind(reason, '''MCsqare''')));
assertFalse(isempty(strfind(reason, 'unknown')));
assertFalse(isempty(strfind(reason, 'SVDPB')));

function test_resolveUnsupportedRadiationMode
% the radiation mode is checked before the availability check of the class,
% which may report misleading messages for other radiation modes
pln = helper_photonPln();
[availableClasses, allClasses] = helper_engineLists(pln);
[ix, reason] = matRad_resolveClassSelection('MCsquare', availableClasses, allClasses, pln);
assertTrue(isempty(ix));
assertFalse(isempty(strfind(reason, '''MCsquare'' does not support radiation mode ''photons''')));
assertFalse(isempty(strfind(reason, 'protons')));

function test_resolveUnsupportedRadiationModeOfStfGenerator
% the availability check of the particle stf generators reports an invalid
% machine for photon machines
pln = helper_photonPln();
availableClasses = matRad_StfGeneratorBase.getAvailableGenerators(pln);
allClasses = matRad_StfGeneratorBase.getAvailableGenerators();
[ix, reason] = matRad_resolveClassSelection('ParticleIMPT', availableClasses, allClasses, pln);
assertTrue(isempty(ix));
assertFalse(isempty(strfind(reason, 'does not support radiation mode ''photons''')));
assertTrue(isempty(strfind(reason, 'invalid')));

function test_resolveUnavailableNameReportsAvailabilityMessage
pln = helper_photonPln();
if DoseEngines.matRad_PhotonVmcEngine.isAvailable(pln)
    moxunit_throw_test_skipped_exception('VMC++ is installed, its availability message cannot be tested');
end
[availableClasses, allClasses] = helper_engineLists(pln);
[ix, reason] = matRad_resolveClassSelection('vmc', availableClasses, allClasses, pln);
assertTrue(isempty(ix));
assertFalse(isempty(strfind(reason, 'VMC++ environment not found')));

function test_resolveListExplainsEveryName
pln = helper_photonPln();
[availableClasses, allClasses] = helper_engineLists(pln);
[ix, reason] = matRad_resolveClassSelection({'MCsquare', 'MCsqare'}, availableClasses, allClasses, pln);
assertTrue(isempty(ix));
assertFalse(isempty(strfind(reason, '''MCsquare''')));
assertFalse(isempty(strfind(reason, '''MCsqare''')));

function test_resolveInvalidSelection
pln = helper_photonPln();
[availableClasses, allClasses] = helper_engineLists(pln);
[ix, reason] = matRad_resolveClassSelection(5, availableClasses, allClasses, pln);
assertTrue(isempty(ix));
assertFalse(isempty(strfind(reason, 'neither a name nor a cell array of names')));
