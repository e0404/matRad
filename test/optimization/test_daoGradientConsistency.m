function test_suite = test_daoGradientConsistency
% Finite-difference consistency checks for the analytic objective gradients
% and constraint Jacobians of the DAO and VMAT optimization problems.
%
% The output should always be test_suite, and the function name the same as
% your file name

% To collect all tests defined below, this is needed in newer Matlab
% versions. test_functions will collect function handles to below test
% functions
test_functions = localfunctions();

% This will initialize the test suite, i.e., take the functions from
% test_functions, check if they contain "test", convert them into a MOxUnit
% Test Case, and add them to the test-runner
initTestSuite;

function cst = helper_prepareCst(cst, dij, pln)
% same cst preparation as matRad_directApertureOptimization: instantiate
% objective/constraint objects, scale to fraction dose, resize to dose grid
cst = matRad_setOverlapPriorities(cst);
for i = 1:size(cst, 1)
    for j = 1:numel(cst{i, 6})
        obj = cst{i, 6}{j};
        if ~isa(obj, 'matRad_DoseOptimizationFunction')
            obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
        end
        obj = obj.setDoseParameters(obj.getDoseParameters() / pln.numOfFractions);
        cst{i, 6}{j} = obj;
    end
end
cst = matRad_resizeCstToGrid(cst, dij.ctGrid.x, dij.ctGrid.y, dij.ctGrid.z, ...
                             dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);

% add one smooth (logsumexp) dose constraint on the target so the
% dosimetric constraint jacobian path is exercised at all (the fixture
% only carries objectives); soften epsilon for finite-difference accuracy
targetIx = find(strcmp(cst(:, 3), 'TARGET'), 1);
dc = DoseConstraints.matRad_MinMaxDose(0, 2.2, 'approx');
dc.epsilon = 0.5;
cst{targetIx, 6}{end + 1} = dc;

function [optiProb, dij, cst, vec] = helper_getDaoProblem()
p = load('photons_testData.mat');
pln = p.pln;
pln.propSeq.sequencer = 'siochi';
resultGUI = matRad_sequencing(p.resultGUI, p.stf, pln);
apertureInfo = resultGUI.sequencing.apertureInfo;
optiProb = matRad_OptimizationProblemDAO(matRad_DoseProjection(), apertureInfo);
dij = p.dij;
cst = helper_prepareCst(p.cst, dij, pln);
vec = helper_perturbVector(apertureInfo);

function [optiProb, dij, cst, vec] = helper_getVmatProblem()
p = load('photons_testData.mat', 'ct', 'cst', 'pln', 'dij');
pln = p.pln;

pln.propStf.generator                = 'PhotonVMAT';
pln.propStf.gantryAngles             = [-180, 180];
pln.propStf.couchAngles              = [0, 0];
pln.propStf.maxGantryAngleSpacing    = 15;
pln.propStf.maxDAOGantryAngleSpacing = 30;
pln.propStf.maxFMOGantryAngleSpacing = 45;
pln.propStf.isoCenter                = matRad_getIsoCenter(p.cst, p.ct, 0);
pln.propSeq.continuousAperture       = false;
pln.propOpt.runVMAT                  = true;

stf = matRad_generateStf(p.ct, p.cst, pln);

sequencer = matRad_SequencingPhotonsSiochiLeaf(pln);
sequencer.runVMAT = true;
sequencer.sequencingLevel = 5;
sequencer.weightToMU = 100;

w = ones(sum([stf.numOfRays]), 1);
sequence = sequencer.sequence(w, stf);
apertureInfo = sequence.apertureInfo;

% synthetic (deterministic) dose influence matrix consistent with the VMAT
% stf: physical plausibility is irrelevant for derivative consistency, but
% the bixels must deposit dose into the cst structures so that the
% dosimetric gradient is nonzero
dij = p.dij;
nVox   = size(p.dij.physicalDose{1}, 1);
nBixel = sum([stf.totalNumOfBixels]);
dij.physicalDose     = {helper_syntheticDose(nVox, nBixel, p.cst)};
dij.totalNumOfBixels = nBixel;
dij.numOfBeams       = numel(stf);

cst = helper_prepareCst(p.cst, dij, pln);
optiProb = matRad_OptimizationProblemVMAT(matRad_DoseProjection(), apertureInfo);
vec = helper_perturbVector(apertureInfo);

function doseMx = helper_syntheticDose(nVox, nBixel, cst)
% deterministic sparse dose matrix: every bixel deposits smoothly varying
% dose into a contiguous block of cst structure voxels
structVox = [];
for s = 1:size(cst, 1)
    structVox = [structVox; cst{s, 4}{1}(:)]; %#ok<AGROW>
end
structVox = unique(structVox);
nStructVox = numel(structVox);

nnzPer = min(40, nStructVox);
rows = zeros(nnzPer * nBixel, 1);
cols = zeros(nnzPer * nBixel, 1);
vals = zeros(nnzPer * nBixel, 1);
idx = 0;
for b = 1:nBixel
    start = mod((b - 1) * 17, nStructVox - nnzPer + 1) + 1;
    rows(idx + (1:nnzPer)) = structVox(start:(start + nnzPer - 1));
    cols(idx + (1:nnzPer)) = b;
    vals(idx + (1:nnzPer)) = 1e-3 * (1.5 + sin(b + (1:nnzPer)'));
    idx = idx + nnzPer;
end
doseMx = sparse(rows, cols, vals, nVox, nBixel);

function vec = helper_perturbVector(apertureInfo)
% Perturb the sequenced aperture vector away from non-differentiable
% points: sequenced leaf positions sit exactly on bixel boundaries (kinks
% of the coverage fraction) and neighboring apertures often share leaf
% positions (kinks of the abs() leaf speed terms). Deterministic,
% index-dependent offsets move every entry to a generic smooth point.
vec = apertureInfo.apertureVector;
nSh = apertureInfo.totalNumOfShapes;
nLP = apertureInfo.totalNumOfLeafPairs;
ixL = nSh + (1:nLP)';
ixR = nSh + nLP + (1:nLP)';

% golden-ratio offsets: distinct for every vector index
off = @(ix) 0.05 + 0.25 * mod(ix * 0.6180339887, 1);

vec(ixL) = vec(ixL) - apertureInfo.bixelWidth * off(ixL);
vec(ixR) = vec(ixR) + apertureInfo.bixelWidth * off(ixR);

% clip into the allowed leaf position range
lo = apertureInfo.limMx(:, 1);
hi = apertureInfo.limMx(:, 2);
vec([ixL; ixR]) = min(max(vec([ixL; ixR]), lo([ixL; ixR]) + 1e-2), hi([ixL; ixR]) - 1e-2);

% strictly positive, pairwise distinct aperture weights
vec(1:nSh) = abs(vec(1:nSh)) .* (1 + 0.1 * off((1:nSh)')) + 1e-2;

% gantry rotation times (VMAT only): keep strictly positive
if numel(vec) > nSh + 2 * nLP
    ixT = ((nSh + 2 * nLP + 1):numel(vec))';
    vec(ixT) = abs(vec(ixT)) .* (1 + 0.1 * off(ixT)) + 1e-2;
end

function sampleIx = helper_sampleIndices(apertureInfo, nVec)
% deterministic sample: a few aperture weights, left/right leaf positions,
% and (if present) time entries
nSh = apertureInfo.totalNumOfShapes;
nLP = apertureInfo.totalNumOfLeafPairs;
sampleIx = [1, max(1, ceil(nSh / 2)), nSh, ...
            nSh + [1, ceil(nLP / 3), nLP], ...
            nSh + nLP + [2, nLP - 1]];
if nVec > nSh + 2 * nLP
    sampleIx = [sampleIx, nSh + 2 * nLP + 1, nVec];
end
sampleIx = unique(sampleIx);

function helper_checkGradient(fFun, g, vec, sampleIx)
anyNonzero = false;
for k = 1:numel(sampleIx)
    i = sampleIx(k);
    h = 1e-4 * max(1, abs(vec(i)));
    vp = vec;
    vp(i) = vp(i) + h;
    vm = vec;
    vm(i) = vm(i) - h;
    gFD = (fFun(vp) - fFun(vm)) / (2 * h);
    tol = 1e-4 * max([1e-2, abs(g(i)), abs(gFD)]);
    assertTrue(abs(g(i) - gFD) <= tol, ...
               sprintf('gradient mismatch at index %d: analytic %g vs FD %g', i, g(i), gFD));
    anyNonzero = anyNonzero || abs(gFD) > 1e-10;
end
% guard against a trivially passing all-zero comparison
assertTrue(anyNonzero);

function helper_checkJacobian(cFun, jacobAn, vec, sampleIx)
anyNonzero = false;
for k = 1:numel(sampleIx)
    i = sampleIx(k);
    h = 1e-4 * max(1, abs(vec(i)));
    vp = vec;
    vp(i) = vp(i) + h;
    vm = vec;
    vm(i) = vm(i) - h;
    colFD = (cFun(vp) - cFun(vm)) / (2 * h);
    colAn = full(jacobAn(:, i));
    tol = 1e-4 * max([1e-2, max(abs(colAn)), max(abs(colFD))]);
    assertTrue(max(abs(colAn - colFD)) <= tol, ...
               sprintf('jacobian mismatch in column %d: max deviation %g', i, max(abs(colAn - colFD))));
    anyNonzero = anyNonzero || any(abs(colFD) > 1e-10);
end
assertTrue(anyNonzero);

function test_daoObjectiveGradientConsistency
[optiProb, dij, cst, vec] = helper_getDaoProblem();
sampleIx = helper_sampleIndices(optiProb.apertureInfo, numel(vec));
fFun = @(v) optiProb.matRad_objectiveFunction(v, dij, cst);
g = optiProb.matRad_objectiveGradient(vec, dij, cst);
helper_checkGradient(fFun, g, vec, sampleIx);

function test_daoConstraintJacobianConsistency
[optiProb, dij, cst, vec] = helper_getDaoProblem();
sampleIx = helper_sampleIndices(optiProb.apertureInfo, numel(vec));
cFun = @(v) optiProb.matRad_constraintFunctions(v, dij, cst);
J = optiProb.matRad_constraintJacobian(vec, dij, cst);
helper_checkJacobian(cFun, J, vec, sampleIx);

function test_vmatObjectiveGradientConsistency
[optiProb, dij, cst, vec] = helper_getVmatProblem();
sampleIx = helper_sampleIndices(optiProb.apertureInfo, numel(vec));
fFun = @(v) optiProb.matRad_objectiveFunction(v, dij, cst);
g = optiProb.matRad_objectiveGradient(vec, dij, cst);
helper_checkGradient(fFun, g, vec, sampleIx);

function test_vmatConstraintJacobianConsistency
[optiProb, dij, cst, vec] = helper_getVmatProblem();
sampleIx = helper_sampleIndices(optiProb.apertureInfo, numel(vec));
cFun = @(v) optiProb.matRad_constraintFunctions(v, dij, cst);
J = optiProb.matRad_constraintJacobian(vec, dij, cst);
helper_checkJacobian(cFun, J, vec, sampleIx);
