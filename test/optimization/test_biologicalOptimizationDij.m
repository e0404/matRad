function test_suite = test_biologicalOptimizationDij

test_functions = localfunctions();

initTestSuite;

function test_variableRBECellAxBxComputesIxDoseAndGamma
[dij, cst] = helper_biologicalDijAndCst(2);

dij = matRad_prepareBiologicalOptimizationDij( ...
                                              dij, cst, matRad_VariableRBEProjection);

assertTrue(iscell(dij.ax));
assertTrue(iscell(dij.bx));
assertTrue(iscell(dij.ixDose));
assertTrue(iscell(dij.gamma));
assertEqual(find(dij.ixDose{1}), [1; 2; 3]);
assertEqual(find(dij.ixDose{2}), [1; 2; 3]);
assertElementsAlmostEqual(dij.gamma{1}(1:3), ...
                          [2; 5 / 3; 2], 'absolute', 1e-12);
assertElementsAlmostEqual(dij.gamma{2}(1:3), ...
                          [5 / 3; 2; 2], 'absolute', 1e-12);

function test_effectCellAxBxComputesIxDoseWithoutGamma
[dij, cst] = helper_biologicalDijAndCst(1);

dij = matRad_prepareBiologicalOptimizationDij( ...
                                              dij, cst, matRad_EffectProjection);

assertTrue(iscell(dij.ixDose));
assertEqual(find(dij.ixDose{1}), [1; 2; 3]);
assertFalse(isfield(dij, 'gamma'));

function test_missingAxBxIsRebuiltFromCst
[dij, cst] = helper_biologicalDijAndCst(2);
dij = rmfield(dij, {'ax', 'bx'});

dij = matRad_prepareBiologicalOptimizationDij( ...
                                              dij, cst, matRad_BEDProjection);

assertTrue(iscell(dij.ax));
assertTrue(iscell(dij.bx));
assertElementsAlmostEqual(dij.ax{1}, [0.2; 0.1; 0.2; 0], ...
                          'absolute', 1e-12);
assertElementsAlmostEqual(dij.bx{2}, [0.03; 0.05; 0.05; 0], ...
                          'absolute', 1e-12);
assertEqual(find(dij.ixDose{1}), [1; 2; 3]);

function test_numericSingleCtAxBxIsNormalizedToCell
[dij, cst] = helper_biologicalDijAndCst(1);
dij.ax = dij.ax{1};
dij.bx = dij.bx{1}';

dij = matRad_prepareBiologicalOptimizationDij( ...
                                              dij, cst, matRad_EffectProjection);

assertTrue(iscell(dij.ax));
assertTrue(iscell(dij.bx));
assertEqual(size(dij.ax{1}), [4 1]);
assertEqual(size(dij.bx{1}), [4 1]);

function test_invalidAxBxSizeFailsClearly
[dij, cst] = helper_biologicalDijAndCst(1);
dij.ax{1} = dij.ax{1}(1:3);

assertExceptionThrown(@() matRad_prepareBiologicalOptimizationDij( ...
                                                                  dij, cst, matRad_EffectProjection), 'matRad:Error');

function [dij, cst] = helper_biologicalDijAndCst(numOfCtScen)
numOfVoxels = 4;
dij = struct();
dij.doseGrid.numOfVoxels = numOfVoxels;
dij.totalNumOfBixels = 2;
dij.numOfScenarios = numOfCtScen;
dij.physicalDose = cell(numOfCtScen, 1, 1);
dij.ax = cell(numOfCtScen, 1);
dij.bx = cell(numOfCtScen, 1);

cst = helper_biologicalCst(numOfCtScen);
for ctScen = 1:numOfCtScen
    doseMatrix = sparse(numOfVoxels, 2);
    doseMatrix(cst{1, 4}{ctScen}, 1) = 1;
    doseMatrix(cst{2, 4}{ctScen}, 2) = 1;
    dij.physicalDose{ctScen, 1, 1} = doseMatrix;
    dij.ax{ctScen} = zeros(numOfVoxels, 1);
    dij.bx{ctScen} = zeros(numOfVoxels, 1);
    dij.ax{ctScen}(cst{1, 4}{ctScen}) = cst{1, 5}.alphaX;
    dij.bx{ctScen}(cst{1, 4}{ctScen}) = cst{1, 5}.betaX;
    dij.ax{ctScen}(cst{2, 4}{ctScen}) = cst{2, 5}.alphaX;
    dij.bx{ctScen}(cst{2, 4}{ctScen}) = cst{2, 5}.betaX;
end

function cst = helper_biologicalCst(numOfCtScen)
cst = cell(2, 6);
cst{1, 1} = 1;
cst{1, 2} = 'TARGET';
cst{1, 3} = 'TARGET';
cst{1, 4} = cell(1, numOfCtScen);
cst{1, 5} = struct('alphaX', 0.2, 'betaX', 0.05);
cst{1, 6} = {};
cst{2, 1} = 2;
cst{2, 2} = 'OAR';
cst{2, 3} = 'OAR';
cst{2, 4} = cell(1, numOfCtScen);
cst{2, 5} = struct('alphaX', 0.1, 'betaX', 0.03);
cst{2, 6} = {};

for ctScen = 1:numOfCtScen
    if ctScen == 1
        cst{1, 4}{ctScen} = [1; 3];
        cst{2, 4}{ctScen} = 2;
    else
        cst{1, 4}{ctScen} = [2; 3];
        cst{2, 4}{ctScen} = 1;
    end
end
