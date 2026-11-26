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

    vheeDummyPln = struct('radiationMode','VHEE','machine','Generic');
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(vheeDummyPln);
    assertTrue(isa(engine,'DoseEngines.matRad_ParticleHongPencilBeamEngine'));

    vheeDummyPlnFocused = struct('radiationMode','VHEE','machine','Focused');
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(vheeDummyPlnFocused);
    assertTrue(isa(engine,'DoseEngines.matRad_ParticleHongPencilBeamEngine'));

function test_getEngineFromPlnByName
protonDummyPln = struct('radiationMode', 'protons', 'machine', 'Generic', 'propDoseCalc', struct('engine', 'MCsquare'));
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(protonDummyPln);
assertTrue(isa(engine, 'DoseEngines.matRad_ParticleMCsquareEngine'));

% Wrong name, falls back to the default engine without strict selection
cleanup = helper_setStrictSelection(false);
photonDummyPln = struct('radiationMode', 'photons', 'machine', 'Generic', 'propDoseCalc', struct('engine', 'MCsquare'));
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(photonDummyPln);
assertTrue(isa(engine, 'DoseEngines.matRad_PhotonPencilBeamSVDEngine'));

function test_getEngineFromPlnReportsUnavailableEngine
% #943: the fallback names the requested engine and why it is not used
cleanup = helper_setStrictSelection(false);
pln = struct('radiationMode', 'photons', 'machine', 'Generic', 'propDoseCalc', struct('engine', 'MCsquare'));
[logEntries, engine] = helper_captureLog(@() DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln));
assertTrue(isa(engine, 'DoseEngines.matRad_PhotonPencilBeamSVDEngine'));
assertTrue(helper_logContains(logEntries, 'WARNING', '''MCsquare'''));
assertTrue(helper_logContains(logEntries, 'WARNING', '''photons'''));
% the fallback engine must not complain about the requested engine
assertFalse(helper_logContains(logEntries, 'WARNING', 'Inconsistent'));

function test_getEngineFromPlnReportsUnknownEngine
cleanup = helper_setStrictSelection(false);
pln = struct('radiationMode', 'photons', 'machine', 'Generic', 'propDoseCalc', struct('engine', 'MCsqare'));
[logEntries, engine] = helper_captureLog(@() DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln));
assertTrue(isa(engine, 'DoseEngines.matRad_PhotonPencilBeamSVDEngine'));
assertTrue(helper_logContains(logEntries, 'WARNING', '''MCsqare'' is unknown'));
assertFalse(helper_logContains(logEntries, 'WARNING', 'Inconsistent'));

function test_getEngineFromPlnReportsAvailabilityMessage
pln = struct('radiationMode', 'photons', 'machine', 'Generic', 'propDoseCalc', struct('engine', 'vmc'));
if DoseEngines.matRad_PhotonVmcEngine.isAvailable(pln)
    moxunit_throw_test_skipped_exception('VMC++ is installed, its availability message cannot be tested');
end
cleanup = helper_setStrictSelection(false);
logEntries = helper_captureLog(@() DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln));
assertTrue(helper_logContains(logEntries, 'WARNING', 'VMC++ environment not found'));

function test_getEngineFromPlnAcceptsListOfEngines
% a list of engines, as in the defaults of MatRad_Config, selects the
% first available engine without a warning
cleanup = helper_setStrictSelection(true);
pln = struct('radiationMode', 'protons', 'machine', 'Generic', 'propDoseCalc', struct('engine', {{'SVDPB', 'HongPB'}}));
[logEntries, engine] = helper_captureLog(@() DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln));
assertTrue(isa(engine, 'DoseEngines.matRad_ParticleHongPencilBeamEngine'));
assertFalse(helper_logContains(logEntries, 'WARNING', 'Inconsistent'));
assertFalse(helper_logContains(logEntries, 'WARNING', 'not valid'));

function test_getEngineFromPlnStrictSelection
cleanup = helper_setStrictSelection(true);

% explicitly requested, but unavailable or unknown engines throw
pln = struct('radiationMode', 'photons', 'machine', 'Generic', 'propDoseCalc', struct('engine', 'MCsquare'));
assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln), 'matRad:Error');
pln.propDoseCalc.engine = 'MCsqare';
assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln), 'matRad:Error');
pln.propDoseCalc.engine = 5;
assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln), 'matRad:Error');

% without a request, the default engine is still used
pln = struct('radiationMode', 'photons', 'machine', 'Generic');
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
assertTrue(isa(engine, 'DoseEngines.matRad_PhotonPencilBeamSVDEngine'));
pln.propDoseCalc = struct('doseGrid', struct('resolution', struct('x', 4, 'y', 4, 'z', 4)));
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
assertTrue(isa(engine, 'DoseEngines.matRad_PhotonPencilBeamSVDEngine'));

function test_engineCaseInsensitiveRequestIsConsistent
% the engine is matched case-insensitively, so this is not inconsistent
cleanup = helper_setStrictSelection(false);
pln = struct('radiationMode', 'photons', 'machine', 'Generic', 'propDoseCalc', struct('engine', 'svdpb'));
[logEntries, engine] = helper_captureLog(@() DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln));
assertTrue(isa(engine, 'DoseEngines.matRad_PhotonPencilBeamSVDEngine'));
assertFalse(helper_logContains(logEntries, 'WARNING', 'Inconsistent'));
