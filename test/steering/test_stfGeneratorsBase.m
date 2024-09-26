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
    assertTrue(isa(generator,'matRad_PhotonStfGeneratorIMRT'));
    assertEqual(generator.radiationMode,'photons');
    assertEqual(generator.machine,'Generic');

    protonDummyPln = struct('radiationMode','protons','machine','Generic');
    generator = matRad_StfGeneratorBase.getGeneratorFromPln(protonDummyPln);
    assertTrue(isa(generator,'matRad_ParticleStfGeneratorIMPT'));
    assertEqual(generator.radiationMode,'protons');
    assertEqual(generator.machine,'Generic');

    carbonDummyPln = struct('radiationMode','carbon','machine','Generic');
    generator = matRad_StfGeneratorBase.getGeneratorFromPln(carbonDummyPln);
    assertTrue(isa(generator,'matRad_ParticleStfGeneratorIMPT'));
    assertEqual(generator.radiationMode,'carbon');
    assertEqual(generator.machine,'Generic');

    heliumDummyPln = struct('radiationMode','helium','machine','Generic');
    generator = matRad_StfGeneratorBase.getGeneratorFromPln(heliumDummyPln);
    assertTrue(isa(generator,'matRad_ParticleStfGeneratorIMPT'));
    assertEqual(generator.radiationMode,'helium');
    assertEqual(generator.machine,'Generic');

    brachyDummyPln = struct('radiationMode','brachy','machine','HDR');
    generator = matRad_StfGeneratorBase.getGeneratorFromPln(brachyDummyPln);
    assertTrue(isa(generator,'matRad_BrachyStfGenerator'));
    assertEqual(generator.radiationMode,'brachy');
    assertEqual(generator.machine,'HDR');

function test_getGeneratorFromPlnByName
    protonDummyPln = struct('radiationMode','protons','machine','Generic','propStf',struct('generator','ParticleIMPT'));
    generator = matRad_StfGeneratorBase.getGeneratorFromPln(protonDummyPln);
    assertTrue(isa(generator,'matRad_ParticleStfGeneratorIMPT'));

    %Wrong name
    photonDummyPln = struct('radiationMode','photons','machine','Generic','propStf',struct('generator','SimpleBrachy'));
    generator = matRad_StfGeneratorBase.getGeneratorFromPln(photonDummyPln);
    assertTrue(isa(generator,'matRad_PhotonStfGeneratorIMRT'));
