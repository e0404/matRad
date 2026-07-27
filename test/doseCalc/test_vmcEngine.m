function test_suite = test_vmcEngine
% Regression tests for the VMC++ DoseEngine integration that do not require
% the separately distributed VMC++ binaries.

test_functions = localfunctions();
initTestSuite;

function test_vmcRootMatchesCurrentSourceLayout
matRad_cfg = MatRad_Config.instance();
expectedRoot = fullfile(matRad_cfg.matRadSrcRoot, 'doseCalc', 'vmc++');

assertEqual(DoseEngines.matRad_PhotonVmcEngine.getVmcRoot(), expectedRoot);

function test_vmcAvailabilityUsesCurrentSourceLayout
pln.radiationMode = 'photons';
pln.machine = 'Generic';
machine = matRad_loadMachine(pln);

[available, ~] = DoseEngines.matRad_PhotonVmcEngine.isAvailable(pln, machine);
expected = isfolder(fullfile(DoseEngines.matRad_PhotonVmcEngine.getVmcRoot(), 'bin'));
assertEqual(available, expected);

function test_vmcRejectsNonCtDoseGrid
ct.cubeDim = [2 3 4];
ct.numOfCtScen = 1;

dij.doseGrid.dimensions = [2 3 3];
dij.doseGrid.numOfVoxels = prod(dij.doseGrid.dimensions);
dij.numOfScenarios = 1;

pln.radiationMode = 'photons';

assertExceptionThrown(@() matRad_calcPhotonDoseVmc(ct, struct(), pln, {}, false, dij), ...
                      'matRad:vmc:InvalidDoseGrid');

function test_vmcEngineRunsCommonInitialization
vmcBin = fullfile(DoseEngines.matRad_PhotonVmcEngine.getVmcRoot(), 'bin');
if isfolder(vmcBin)
    % Do not launch a real Monte Carlo calculation from this unit test.
    return
end

p = load('photons_testData.mat', 'ct', 'cst', 'pln', 'stf');
engine = DoseEngines.matRad_PhotonVmcEngine(p.pln);

assertExceptionThrown(@() engine.calcDoseInfluence(p.ct, p.cst, p.stf), ...
                      'matRad:vmc:EnvironmentNotFound');
