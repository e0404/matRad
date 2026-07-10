function test_suite = test_sigmaRashi

test_functions = localfunctions();

initTestSuite;

function test_calcSigmaRashi

baseDataEntry.range = 100; % [mm]

rangeShifter.ID = 1;
rangeShifter.eqThickness = 1;
rangeShifter.sourceRashiDistance = 9000;

SSD = 10000;

sigma = matRad_calcSigmaRashi(baseDataEntry, rangeShifter, SSD);

assertElementsAlmostEqual(sigma, 4.023, 'relative', 1e-2, 1e-2);

% explicit radiation mode gives the same result
sigmaExplicit = matRad_calcSigmaRashi(baseDataEntry, rangeShifter, SSD, 'protons');
assertElementsAlmostEqual(sigmaExplicit, sigma);

function test_calcSigmaRashiFromEnergy
% energy-based fallback must be consistent with the range-based
% computation (Bragg-Kleeman with energy per nucleon)
alpha = 0.0022;
p = 1.77;

rangeShifter.ID = 1;
rangeShifter.eqThickness = 1;
rangeShifter.sourceRashiDistance = 9000;

SSD = 10000;

baseDataEntry.range = 100; % [mm]

% proton with 10 cm range
energyEntry.energy = (10 / alpha)^(1 / p);
assertElementsAlmostEqual(matRad_calcSigmaRashi(energyEntry, rangeShifter, SSD), ...
                          matRad_calcSigmaRashi(baseDataEntry, rangeShifter, SSD), 'relative', 1e-6);

% carbon with 10 cm range (equivalent proton range is z^2/A = 3 times larger)
energyEntry.energy = (3 * 10 / alpha)^(1 / p);
assertElementsAlmostEqual(matRad_calcSigmaRashi(energyEntry, rangeShifter, SSD, 'carbon'), ...
                          matRad_calcSigmaRashi(baseDataEntry, rangeShifter, SSD, 'carbon'), 'relative', 1e-6);

function test_calcSigmaRashiIons

baseDataEntry.range = 100; % [mm]

rangeShifter.ID = 1;
rangeShifter.eqThickness = 1;
rangeShifter.sourceRashiDistance = 9000;

SSD = 10000;

sigmaP  = matRad_calcSigmaRashi(baseDataEntry, rangeShifter, SSD, 'protons');
sigmaHe = matRad_calcSigmaRashi(baseDataEntry, rangeShifter, SSD, 'helium');
sigmaC  = matRad_calcSigmaRashi(baseDataEntry, rangeShifter, SSD, 'carbon');
sigmaO  = matRad_calcSigmaRashi(baseDataEntry, rangeShifter, SSD, 'oxygen');

% helium (z^2/A = 1) scatters exactly half as much as a proton with
% the same range
assertElementsAlmostEqual(sigmaHe, 0.5 * sigmaP, 'relative', 1e-10);

% heavier ions scatter less at the same range
assertTrue(sigmaP > sigmaHe && sigmaHe > sigmaC && sigmaC > sigmaO);

% carbon-to-proton scattering ratio at same range is about 0.27
assertElementsAlmostEqual(sigmaC / sigmaP, 0.269, 'relative', 1e-2);

function test_calcSigmaRashiInvalidInput

baseDataEntry.range = 100; % [mm]

rangeShifter.ID = 1;
rangeShifter.eqThickness = 1;
rangeShifter.sourceRashiDistance = 9000;

SSD = 10000;

% unsupported radiation mode
assertExceptionThrown(@() matRad_calcSigmaRashi(baseDataEntry, rangeShifter, SSD, 'VHEE'));

% range shifter too thick
thickRashi = rangeShifter;
thickRashi.eqThickness = 99;
assertExceptionThrown(@() matRad_calcSigmaRashi(baseDataEntry, thickRashi, SSD));
