function test_suite = test_photonCollapsedConeEngine

test_functions = localfunctions();

initTestSuite;

function [ct, cst, pln, stf] = helper_waterBox(numVox, res, targetSize, bixelWidth)
builder = matRad_PhantomBuilder([numVox numVox numVox], [res res res], 1);
builder.addBoxTarget('Target', targetSize * [1 1 1], 'HU', 0, ...
                     'objectives', struct(DoseObjectives.matRad_SquaredDeviation(800, 60)));
builder.addBoxOAR('Water', (numVox - 4) * [1 1 1], 'HU', 0, ...
                  'objectives', struct(DoseObjectives.matRad_SquaredOverdosing(10, 30)));
[ct, cst] = builder.getctcst();

pln.radiationMode = 'photons';
pln.machine = 'Generic';
pln.numOfFractions = 1;
pln.propStf.gantryAngles = 0;
pln.propStf.couchAngles = 0;
pln.propStf.bixelWidth = bixelWidth;
pln.propStf.isoCenter = matRad_getIsoCenter(cst, ct, 0);
pln.propDoseCalc.engine = 'CCC';
pln.propDoseCalc.doseGrid.resolution = struct('x', res, 'y', res, 'z', res);
pln.propDoseCalc.numPolarAngles = 4;
pln.propDoseCalc.numAzimuthalAngles = 4;
pln.propDoseCalc.scatterCutOff = 50;
pln.propDoseCalc.enableDijSampling = false;
stf = matRad_generateStf(ct, cst, pln);

function fileName = helper_shortRangeKernel(numPolarAngles)
edges = linspace(0, 180, numPolarAngles + 1);
cosEdges = cosd(edges);
solidAngle = (cosEdges(1:end - 1) - cosEdges(2:end)) / 2;
ccKernel.fit.numPolarAngles = numPolarAngles;
ccKernel.fit.polarEdges = edges;
ccKernel.fit.polarAngles = acosd(0.5 * (cosEdges(1:end - 1) + cosEdges(2:end)));
ccKernel.fit.A = 0.6 * solidAngle;
ccKernel.fit.a = 5 * ones(1, numPolarAngles);
ccKernel.fit.B = 0.3 * solidAngle;
ccKernel.fit.b = 3 * ones(1, numPolarAngles);
fileName = fullfile(helper_temporaryFolder('matRad_test_CCC'), 'shortRangeKernel.mat');
save(fileName, 'ccKernel');

function test_getCCCEngineFromPln
pln = struct('radiationMode', 'photons', 'machine', 'Generic');
pln.propDoseCalc.engine = 'CCC';
pln.propDoseCalc.numAzimuthalAngles = 6;
pln.propDoseCalc.scatterCutOff = 80;
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
assertTrue(isa(engine, 'DoseEngines.matRad_PhotonCollapsedConeEngine'));
assertEqual(engine.numAzimuthalAngles, 6);
assertEqual(engine.scatterCutOff, 80);
assertEqual(engine.numPolarAngles, 12);
assertEqual(engine.kernelTilt, 'beam');

function test_isAvailable
pln = struct('radiationMode', 'photons', 'machine', 'Generic');
assertTrue(DoseEngines.matRad_PhotonCollapsedConeEngine.isAvailable(pln));
pln = struct('radiationMode', 'protons', 'machine', 'Generic');
assertFalse(DoseEngines.matRad_PhotonCollapsedConeEngine.isAvailable(pln));
pln = struct('radiationMode', 'photons', 'machine', 'Generic');
pln.propDoseCalc.kernelFile = fullfile(tempdir(), 'matRad_noSuchKernel.mat');
assertFalse(DoseEngines.matRad_PhotonCollapsedConeEngine.isAvailable(pln));

function test_kernelData
fileName = DoseEngines.matRad_PhotonCollapsedConeEngine.defaultKernelFile();
kernel = DoseEngines.matRad_PhotonCollapsedConeEngine.loadKernel(fileName, 12);
assertEqual(numel(kernel.A), 12);
assertTrue(all(kernel.A >= 0 & kernel.B >= 0 & kernel.a > 0 & kernel.b > 0));
total = sum(kernel.A + kernel.B);
assertTrue(total > 0.9 && total < 1);
kernel6 = DoseEngines.matRad_PhotonCollapsedConeEngine.loadKernel(fileName, 6);
assertElementsAlmostEqual(sum(kernel6.A + kernel6.B), total, 'relative', 1e-10);
kernel2 = DoseEngines.matRad_PhotonCollapsedConeEngine.loadKernel(fileName, 12, 2);
assertEqual(size(kernel2.B), [2 12]);
assertElementsAlmostEqual(sum(kernel2.A) + sum(kernel2.B(:)), total, 'relative', 1e-10);

function test_dijSmallBox
[ct, cst, pln, stf] = helper_waterBox(16, 6, 4, 12);
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
assertEqual(size(dij.physicalDose{1}), [dij.doseGrid.numOfVoxels sum([stf.totalNumOfBixels])]);
D = dij.physicalDose{1};
assertTrue(nnz(D) > 0);
assertTrue(all(isfinite(nonzeros(D))));
assertTrue(all(nonzeros(D) >= 0));
assertTrue(all(full(sum(D, 1)) > 0));

function test_forwardEqualsDijTimesWeights
[ct, cst, pln, stf] = helper_waterBox(16, 6, 4, 12);
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
w = linspace(0.5, 1.5, dij.totalNumOfBixels)';
resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, w);
assertEqual(dij.doseGrid.dimensions, dij.ctGrid.dimensions);
doseDij = reshape(full(dij.physicalDose{1} * w), dij.doseGrid.dimensions);
assertElementsAlmostEqual(resultGUI.physicalDose, doseDij, 'absolute', 1e-10 * max(doseDij(:)));
assertTrue(all(resultGUI.physicalDose(:) >= 0));
assertFalse(any(isnan(resultGUI.physicalDose(:))));

function test_mexMatchesMatlabEngine
if exist('matRad_collapsedConeMex', 'file') ~= 3
    moxunit_throw_test_skipped_exception('matRad_collapsedConeMex not compiled');
end
[ct, cst, pln, stf] = helper_waterBox(16, 6, 4, 12);
ct.cubeHU{1}(:, 7:9, :) = -700;
w = ones(sum([stf.totalNumOfBixels]), 1);
pln.propDoseCalc.useMex = true;
doseMex = matRad_calcDoseForward(ct, cst, stf, pln, w);
pln.propDoseCalc.useMex = false;
doseRef = matRad_calcDoseForward(ct, cst, stf, pln, w);
relTol = 1e-10;
assertElementsAlmostEqual(doseMex.physicalDose, doseRef.physicalDose, 'absolute', relTol * max(doseRef.physicalDose(:)));

function test_shortRangeKernelCentralAxis
[ct, cst, pln, stf] = helper_waterBox(30, 6, 10, 20);
pln.propDoseCalc.kernelFile = helper_shortRangeKernel(4);
pln.propDoseCalc.scatterCutOff = 20;
engine = DoseEngines.matRad_PhotonCollapsedConeEngine(pln);
resultGUI = engine.calcDoseForward(ct, cst, stf, ones(sum([stf.totalNumOfBixels]), 1));
assertTrue(all(resultGUI.physicalDose(:) >= 0));
assertFalse(any(isnan(resultGUI.physicalDose(:))));

machine = matRad_loadMachine(pln);
center = round(ct.cubeDim / 2);
column = resultGUI.physicalDose(:, center(2), center(3));
first = find(column > 0, 1);
rows = (first + 3:first + 20)';
profile = column(rows);
depth = (rows - first + 0.5) * ct.resolution.y;
voxIx = sub2ind(ct.cubeDim, rows, center(2) * ones(numel(rows), 1), center(3) * ones(numel(rows), 1));
voxWorld = matRad_cubeIndex2worldCoords(voxIx, ct);
source = stf.isoCenter + stf.sourcePoint;
invSquare = (machine.meta.SAD ./ sqrt(sum((voxWorld - repmat(source, numel(rows), 1)).^2, 2))).^2;

terma = DoseEngines.matRad_PhotonCollapsedConeEngine.termaDepthFactor(engine.hardeningMu0, engine.hardeningBeta, depth) .* invSquare;
expected = engine.normalizationFactor * sum(engine.kernel.A + engine.kernel.B) * terma;
assertElementsAlmostEqual(profile, expected, 'relative', 0.01);

dose = interp1(depth, profile ./ invSquare, engine.normalizationDepth);
reference = DoseEngines.matRad_PhotonCollapsedConeEngine.calcAnalyticalDepthDose(machine, engine.normalizationDepth);
assertElementsAlmostEqual(dose, reference, 'relative', 0.01);

function test_hardeningFit
machine = matRad_loadMachine(struct('radiationMode', 'photons', 'machine', 'Generic'));
fileName = DoseEngines.matRad_PhotonCollapsedConeEngine.defaultKernelFile();
kernel = DoseEngines.matRad_PhotonCollapsedConeEngine.loadKernel(fileName, 12);
[mu0, beta] = DoseEngines.matRad_PhotonCollapsedConeEngine.fitHardening(kernel, machine, 5, 100);
assertTrue(mu0 > 0 && beta > 0);
depths = [30 100 200 300];
reference = DoseEngines.matRad_PhotonCollapsedConeEngine.calcAnalyticalDepthDose(machine, depths);
reference = reference / reference(2);
corrected = DoseEngines.matRad_PhotonCollapsedConeEngine.calcBroadFieldDepthDose(kernel, mu0, depths, 5, beta);
corrected = corrected / corrected(2);
mono = DoseEngines.matRad_PhotonCollapsedConeEngine.calcBroadFieldDepthDose(kernel, machine.data.m, depths, 5, 0);
mono = mono / mono(2);
assertTrue(max(abs(corrected ./ reference - 1)) < 0.02);
assertTrue(max(abs(corrected ./ reference - 1)) < max(abs(mono ./ reference - 1)));
assertElementsAlmostEqual(DoseEngines.matRad_PhotonCollapsedConeEngine.termaDepthFactor(0.005, 0, depths), ...
                          0.005 * exp(-0.005 * depths), 'relative', 1e-14);

function test_hardeningCorrectionDisabled
[ct, cst, pln, stf] = helper_waterBox(12, 6, 4, 12);
pln.propDoseCalc.hardeningCorrection = false;
engine = DoseEngines.matRad_PhotonCollapsedConeEngine(pln);
engine.calcDoseInfluence(ct, cst, stf);
machine = matRad_loadMachine(pln);
assertEqual(engine.hardeningBeta, 0);
assertEqual(engine.hardeningMu0, machine.data.m);

function test_kernelEnergyMismatch
pln = struct('radiationMode', 'photons', 'machine', 'Generic');
machine = matRad_loadMachine(pln);
fileName = DoseEngines.matRad_PhotonCollapsedConeEngine.defaultKernelFile();
data = load(fileName, 'ccKernel');
[ok, msg] = DoseEngines.matRad_PhotonCollapsedConeEngine.checkKernelEnergy(data.ccKernel.meta, machine);
assertTrue(ok);
assertTrue(isempty(msg));

machine.data.energy = 10;
assertFalse(DoseEngines.matRad_PhotonCollapsedConeEngine.isAvailable(pln, machine));

ccKernel = data.ccKernel;
ccKernel.meta.nominalEnergy = 10;
kernelFile = fullfile(helper_temporaryFolder('matRad_test_CCC'), 'kernel10MV.mat');
save(kernelFile, 'ccKernel');
[ct, cst, pln, stf] = helper_waterBox(12, 6, 4, 12);
pln.propDoseCalc.kernelFile = kernelFile;
engine = DoseEngines.matRad_PhotonCollapsedConeEngine(pln);
assertExceptionThrown(@() engine.calcDoseInfluence(ct, cst, stf));

function test_invalidKernelTilt
[ct, cst, pln, stf] = helper_waterBox(12, 6, 4, 12);
pln.propDoseCalc.kernelTilt = 'voxel';
assertExceptionThrown(@() matRad_calcDoseInfluence(ct, cst, stf, pln));
