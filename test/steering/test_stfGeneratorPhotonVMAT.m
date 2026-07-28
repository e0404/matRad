function test_suite = test_stfGeneratorPhotonVMAT
% Geometry tests for the VMAT stf generator's arc subdivision.
%
% The generator partitions each arc into nested FMO / DAO / dose sectors.
% The max*GantryAngleSpacing properties are upper bounds, so the derivation
% may go finer than requested but must never go coarser, and must never trade
% away apertures per fluence map to satisfy them.

test_functions = localfunctions();

initTestSuite;

function pln = helper_vmatPln(maxG, maxDAO, maxFMO, gantryAngles, continuousAperture)
if nargin < 4 || isempty(gantryAngles)
    gantryAngles = [-180, 180];
end
if nargin < 5
    continuousAperture = false;
end
p = load('photons_testData.mat', 'ct', 'cst', 'pln');
pln = p.pln;
pln.propStf.generator                = 'PhotonVMAT';
pln.propStf.gantryAngles             = gantryAngles;
pln.propStf.couchAngles              = zeros(1, numel(gantryAngles));
pln.propStf.maxGantryAngleSpacing    = maxG;
pln.propStf.maxDAOGantryAngleSpacing = maxDAO;
pln.propStf.maxFMOGantryAngleSpacing = maxFMO;
pln.propStf.isoCenter                = matRad_getIsoCenter(p.cst, p.ct, 0);
pln.propSeq.continuousAperture       = continuousAperture;
pln.propOpt.runVMAT                  = true;

function stf = helper_generate(pln)
p = load('photons_testData.mat', 'ct', 'cst');
stf = matRad_generateStf(p.ct, p.cst, pln);

function helper_assertGeometryInvariants(stf, maxG, maxDAO, maxFMO, minApert)
gantryAngles = [stf.gantryAngle];
isDAO = arrayfun(@(s) s.arc.isDAOBeam, stf);
isFMO = arrayfun(@(s) s.arc.isFMOBeam, stf);

% nesting: FMO angles are DAO angles are dose angles
assertTrue(all(ismember(find(isFMO), find(isDAO))), 'FMO beams must be DAO beams');
assertTrue(nnz(isFMO) >= 1 && nnz(isDAO) >= 1);

% no two beams on the same physical gantry position (closed-arc endpoint)
assertEqual(numel(gantryAngles), numel(unique(round(mod(gantryAngles, 360), 6))));

% realized spacings never exceed the requested maxima
daoAngles = gantryAngles(isDAO);
fmoAngles = gantryAngles(isFMO);
assertTrue(abs(gantryAngles(2) - gantryAngles(1)) <= maxG + 1e-9);
assertTrue(abs(daoAngles(2) - daoAngles(1)) <= maxDAO + 1e-9);
assertTrue(abs(fmoAngles(2) - fmoAngles(1)) <= maxFMO + 1e-9);

% uniform spacing at every level
assertElementsAlmostEqual(diff(gantryAngles), ...
                          repmat(gantryAngles(2) - gantryAngles(1), 1, numel(gantryAngles) - 1));
assertElementsAlmostEqual(diff(daoAngles), repmat(daoAngles(2) - daoAngles(1), 1, numel(daoAngles) - 1));

% every fluence map gets the same, odd, sufficient aperture budget
kids = arrayfun(@(s) s.arc.numOfChildren, stf(isFMO));
assertTrue(all(kids == kids(1)), 'aperture budget must be identical for every FMO beam');
assertTrue(kids(1) >= minApert, sprintf('only %d apertures per fluence map', kids(1)));
assertEqual(mod(kids(1), 2), 1, 'aperture budget must be odd so the FMO angle stays centred');
assertEqual(nnz(isDAO), nnz(isFMO) * kids(1));

% dose sector bookkeeping is finite and positive everywhere
for i = 1:numel(stf)
    assertTrue(all(isfinite(stf(i).arc.doseAngleBorders)));
    assertTrue(stf(i).arc.doseAngleBordersDiff > 0);
end

function test_vmatGeometryInvariantsAcrossConfigurations
cfgs = [15 30 45; 15 30 90; 4 8 32; 5 10 70; 30 60 180];
for r = 1:size(cfgs, 1)
    for continuousAperture = [false true]
        pln = helper_vmatPln(cfgs(r, 1), cfgs(r, 2), cfgs(r, 3), [], continuousAperture);
        stf = helper_generate(pln);
        helper_assertGeometryInvariants(stf, cfgs(r, 1), cfgs(r, 2), cfgs(r, 3), 3);
    end
end

function test_vmatGeometryInvariantsForReverseAndPartialArcs
arcs = {[180, -180], [0, 180], [30, -150]};
for a = 1:numel(arcs)
    for continuousAperture = [false true]
        pln = helper_vmatPln(15, 30, 45, arcs{a}, continuousAperture);
        stf = helper_generate(pln);
        helper_assertGeometryInvariants(stf, 15, 30, 45, 3);
        expectedSign = sign(arcs{a}(end) - arcs{a}(1));
        assertTrue(all(sign(diff([stf.gantryAngle])) == expectedSign));
    end
end

function test_vmatFormerlyDegenerateConfigYieldsModulation
% Regression: with the previous floor-based derivation, 15/30/45 collapsed
% to a single aperture per fluence map, which strips all modulation out of
% the sequenced plan and leaves DAO nothing to recover. Any maxFMO from the
% realized DAO spacing up to 3x it behaved the same way, so this pins the
% config that example 22 shipped with when the collapse was found.
pln = helper_vmatPln(15, 30, 45);
stf = helper_generate(pln);

isFMO = arrayfun(@(s) s.arc.isFMOBeam, stf);
kids = arrayfun(@(s) s.arc.numOfChildren, stf(isFMO));
assertTrue(all(kids >= 3), sprintf('only %d apertures per fluence map', min(kids)));

% pin the realized geometry so it cannot drift unnoticed
assertEqual(numel(stf), 24);
assertEqual(nnz(arrayfun(@(s) s.arc.isDAOBeam, stf)), 24);
assertEqual(nnz(isFMO), 8);
assertEqual(kids(1), 3);
assertElementsAlmostEqual([stf.gantryAngle], -172.5:15:172.5);

function test_vmatBeamsStayStrictlyInsideTheArc
% Centre sampling: every beam sits at the centre of its dose sector, so no
% beam lands on an arc boundary.
pln = helper_vmatPln(15, 30, 45, [0, 180]);
stf = helper_generate(pln);
gantryAngles = [stf.gantryAngle];

assertTrue(all(gantryAngles > 0) && all(gantryAngles < 180));
halfSector = (gantryAngles(2) - gantryAngles(1)) / 2;
assertElementsAlmostEqual(gantryAngles(1), 0 + halfSector);
assertElementsAlmostEqual(gantryAngles(end), 180 - halfSector);

function test_vmatMinAperturesPerFMOBeamIsHonoured
pln = helper_vmatPln(15, 30, 45);
pln.propStf.minAperturesPerFMOBeam = 5;
stf = helper_generate(pln);

isFMO = arrayfun(@(s) s.arc.isFMOBeam, stf);
kids = arrayfun(@(s) s.arc.numOfChildren, stf(isFMO));
assertTrue(all(kids >= 5));
helper_assertGeometryInvariants(stf, 15, 30, 45, 5);

function test_vmatExplicitSubdivisionFactorsOverrideTheBounds
pln = helper_vmatPln(15, 30, 45);
pln.propStf.aperturesPerFMOBeam  = 5;
pln.propStf.doseAnglesPerDAOBeam = 3;
stf = helper_generate(pln);

isFMO = arrayfun(@(s) s.arc.isFMOBeam, stf);
isDAO = arrayfun(@(s) s.arc.isDAOBeam, stf);
kids = arrayfun(@(s) s.arc.numOfChildren, stf(isFMO));

assertEqual(kids(1), 5);
assertEqual(nnz(isDAO), nnz(isFMO) * 5);
assertEqual(numel(stf), nnz(isDAO) * 3);

function test_vmatRejectsEvenSubdivisionFactors
pln = helper_vmatPln(15, 30, 45);
pln.propStf.aperturesPerFMOBeam = 4;
assertExceptionThrown(@() helper_generate(pln));

pln = helper_vmatPln(15, 30, 45);
pln.propStf.doseAnglesPerDAOBeam = 2;
assertExceptionThrown(@() helper_generate(pln));

function test_vmatRejectsRunawayDAOAngleCount
% Because the spacings are upper bounds, a tight FMO bound with a loose DAO
% bound legitimately forces a very fine DAO grid - that must fail loudly
% rather than silently build an enormous optimization problem.
pln = helper_vmatPln(15, 30, 45);
pln.propStf.maxNumOfDAOAngles = 10;
assertExceptionThrown(@() helper_generate(pln));

function test_vmatRejectsInconsistentSpacingBounds
pln = helper_vmatPln(30, 15, 45);   % maxG > maxDAO
assertExceptionThrown(@() helper_generate(pln));
