function test_suite = test_stfGeneratorsBase

test_functions=localfunctions();

initTestSuite;

function test_stfGeneratorBaseIsAbstract
    if moxunit_util_platform_is_octave()
        assertExceptionThrown(@() matRad_StfGeneratorBase(),'');
    else
        assertExceptionThrown(@() matRad_StfGeneratorBase(),'MATLAB:class:abstract');
    end

function test_abstractIsAvailable
    assertExceptionThrown(@() matRad_StfGeneratorBase.isAvailable,'matRad:Error');

function test_getAvailableGenerators
    avail = matRad_StfGeneratorBase.getAvailableGenerators();
    assertTrue(~isempty(avail));
    assertTrue(isstruct(avail));
    assertTrue(iscolumn(avail));

    photonDummyPln = struct('radiationMode','photons','machine','Generic');
    availPhotons = matRad_StfGeneratorBase.getAvailableGenerators(photonDummyPln);
    assertTrue(~isempty(avail));
    assertTrue(isstruct(avail));
    assertTrue(iscolumn(avail));
    assertTrue(numel(availPhotons) < numel(avail));

function test_loadMachine
    machine = matRad_StfGeneratorBase.loadMachine('photons','Generic');
    assertTrue(isstruct(machine));

    if moxunit_util_platform_is_octave()
        assertExceptionThrown(@() matRad_StfGeneratorBase.loadMachine());
        assertExceptionThrown(@() matRad_StfGeneratorBase.loadMachine('photons'));
    else
        assertExceptionThrown(@() matRad_StfGeneratorBase.loadMachine(),'MATLAB:minrhs');
        assertExceptionThrown(@() matRad_StfGeneratorBase.loadMachine('photons'),'MATLAB:minrhs');
    end
    assertExceptionThrown(@() matRad_StfGeneratorBase.loadMachine('grbl','grbl'),'matRad:Error');

function test_getGeneratorFromPlnDefaults
    photonDummyPln = struct('radiationMode','photons','machine','Generic');
    generator = matRad_StfGeneratorBase.getGeneratorFromPln(photonDummyPln);
    assertTrue(isa(generator,'matRad_StfGeneratorPhotonIMRT'));
    assertEqual(generator.radiationMode,'photons');
    assertEqual(generator.machine,'Generic');

    protonDummyPln = struct('radiationMode','protons','machine','Generic');
    generator = matRad_StfGeneratorBase.getGeneratorFromPln(protonDummyPln);
    assertTrue(isa(generator,'matRad_StfGeneratorParticleIMPT'));
    assertEqual(generator.radiationMode,'protons');
    assertEqual(generator.machine,'Generic');

    carbonDummyPln = struct('radiationMode','carbon','machine','Generic');
    generator = matRad_StfGeneratorBase.getGeneratorFromPln(carbonDummyPln);
    assertTrue(isa(generator,'matRad_StfGeneratorParticleIMPT'));
    assertEqual(generator.radiationMode,'carbon');
    assertEqual(generator.machine,'Generic');

    heliumDummyPln = struct('radiationMode','helium','machine','Generic');
    generator = matRad_StfGeneratorBase.getGeneratorFromPln(heliumDummyPln);
    assertTrue(isa(generator,'matRad_StfGeneratorParticleIMPT'));
    assertEqual(generator.radiationMode,'helium');
    assertEqual(generator.machine,'Generic');

    brachyDummyPln = struct('radiationMode','brachy','machine','HDR');
    generator = matRad_StfGeneratorBase.getGeneratorFromPln(brachyDummyPln);
    assertTrue(isa(generator,'matRad_StfGeneratorBrachy'));
    assertEqual(generator.radiationMode,'brachy');
    assertEqual(generator.machine,'HDR');

function test_getGeneratorFromPlnByName
protonDummyPln = struct('radiationMode', 'protons', 'machine', 'Generic', 'propStf', struct('generator', 'ParticleIMPT'));
generator = matRad_StfGeneratorBase.getGeneratorFromPln(protonDummyPln);
assertTrue(isa(generator, 'matRad_StfGeneratorParticleIMPT'));

% Wrong name, falls back to the default generator without strict selection
cleanup = helper_setStrictSelection(false);
photonDummyPln = struct('radiationMode', 'photons', 'machine', 'Generic', 'propStf', struct('generator', 'SimpleBrachy'));
generator = matRad_StfGeneratorBase.getGeneratorFromPln(photonDummyPln);
assertTrue(isa(generator, 'matRad_StfGeneratorPhotonIMRT'));

function test_getGeneratorFromPlnReportsUnavailableGenerator
% #943: the fallback names the requested generator and why it is not used
cleanup = helper_setStrictSelection(false);
pln = struct('radiationMode', 'photons', 'machine', 'Generic', 'propStf', struct('generator', 'ParticleIMPT'));
[logEntries, generator] = helper_captureLog(@() matRad_StfGeneratorBase.getGeneratorFromPln(pln));
assertTrue(isa(generator, 'matRad_StfGeneratorPhotonIMRT'));
assertTrue(helper_logContains(logEntries, 'WARNING', '''ParticleIMPT'''));
assertTrue(helper_logContains(logEntries, 'WARNING', '''photons'''));
assertFalse(helper_logContains(logEntries, 'WARNING', 'Inconsistent'));

function test_getGeneratorFromPlnReportsUnknownGeneratorWithoutDefaultWarning
% warnDefault only silences the use of the default, not the reason
cleanup = helper_setStrictSelection(false);
pln = struct('radiationMode', 'photons', 'machine', 'Generic', 'propStf', struct('generator', 'PhotonIMRTT'));
[logEntries, generator] = helper_captureLog(@() matRad_StfGeneratorBase.getGeneratorFromPln(pln, false));
assertTrue(isa(generator, 'matRad_StfGeneratorPhotonIMRT'));
assertTrue(helper_logContains(logEntries, 'WARNING', '''PhotonIMRTT'' is unknown'));

function test_getGeneratorFromPlnStrictSelection
cleanup = helper_setStrictSelection(true);
pln = struct('radiationMode', 'photons', 'machine', 'Generic', 'propStf', struct('generator', 'ParticleIMPT'));
assertExceptionThrown(@() matRad_StfGeneratorBase.getGeneratorFromPln(pln), 'matRad:Error');
pln.propStf.generator = 'PhotonIMRTT';
assertExceptionThrown(@() matRad_StfGeneratorBase.getGeneratorFromPln(pln), 'matRad:Error');

% a list of generators, as in the defaults, and no request still work
pln.propStf.generator = {'PhotonIMRT', 'ParticleIMPT'};
generator = matRad_StfGeneratorBase.getGeneratorFromPln(pln);
assertTrue(isa(generator, 'matRad_StfGeneratorPhotonIMRT'));
pln = struct('radiationMode', 'photons', 'machine', 'Generic');
generator = matRad_StfGeneratorBase.getGeneratorFromPln(pln);
assertTrue(isa(generator, 'matRad_StfGeneratorPhotonIMRT'));
