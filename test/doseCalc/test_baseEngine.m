function test_suite = test_baseEngine

test_functions=localfunctions();

initTestSuite;

function test_doseEngineBaseAbstract
    if moxunit_util_platform_is_octave()
        assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase(),'');
    else
        assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase(),'MATLAB:class:abstract');
    end

function test_abstractIsAvailable
    assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase.isAvailable,'matRad:Error');

function test_getAvailableEngines
    avail = DoseEngines.matRad_DoseEngineBase.getAvailableEngines();
    assertTrue(~isempty(avail));
    assertTrue(isstruct(avail));
    assertTrue(iscolumn(avail));

    photonDummyPln = struct('radiationMode','photons','machine','Generic');
    availPhotons = DoseEngines.matRad_DoseEngineBase.getAvailableEngines(photonDummyPln);
    assertTrue(~isempty(avail));
    assertTrue(isstruct(avail));
    assertTrue(iscolumn(avail));
    assertTrue(numel(availPhotons) < numel(avail));

function test_loadMachine
    machine = DoseEngines.matRad_DoseEngineBase.loadMachine('photons','Generic');
    assertTrue(isstruct(machine));

    if moxunit_util_platform_is_octave()
        assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase.loadMachine());
        assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase.loadMachine('photons'));
    else
        assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase.loadMachine(),'MATLAB:minrhs');
        assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase.loadMachine('photons'),'MATLAB:minrhs');
    end
    assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase.loadMachine('grbl','grbl'),'matRad:Error');

function test_getEngineFromPlnDefaults
    photonDummyPln = struct('radiationMode','photons','machine','Generic');
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(photonDummyPln);
    assertTrue(isa(engine,'DoseEngines.matRad_PhotonPencilBeamSVDEngine'));

    protonDummyPln = struct('radiationMode','protons','machine','Generic');
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(protonDummyPln);
    assertTrue(isa(engine,'DoseEngines.matRad_ParticleHongPencilBeamEngine'));

    carbonDummyPln = struct('radiationMode','carbon','machine','Generic');
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(carbonDummyPln);
    assertTrue(isa(engine,'DoseEngines.matRad_ParticleHongPencilBeamEngine'));

    heliumDummyPln = struct('radiationMode','helium','machine','Generic');
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(heliumDummyPln);
    assertTrue(isa(engine,'DoseEngines.matRad_ParticleHongPencilBeamEngine'));

function test_getEngineFromPlnByName
    protonDummyPln = struct('radiationMode','protons','machine','Generic','propDoseCalc',struct('engine','MCsquare'));
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(protonDummyPln);
    assertTrue(isa(engine,'DoseEngines.matRad_ParticleMCsquareEngine'));

    %Wrong name
    photonDummyPln = struct('radiationMode','photons','machine','Generic','propDoseCalc',struct('engine','MCsquare'));
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(photonDummyPln);
    assertTrue(isa(engine,'DoseEngines.matRad_PhotonPencilBeamSVDEngine'));

    


