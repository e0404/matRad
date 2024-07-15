function test_suite = test_BEDprojection

test_functions = localfunctions();

initTestSuite;

function test_BEDprojectionConstructor
bed = matRad_BEDProjection;
assertTrue(isobject(bed));
assertTrue(isa(bed, 'matRad_BEDProjection'));

%test compute single scenario BED projection 
function test_BED_computeSingleScenario
dij.physicalDose{1} = 2;
dij.mAlphaDose{1} = 1;
dij.mSqrtBetaDose{1} = sqrt(0.05)*2;
dij.ax{1} = 0.5;
dij.bx{1} = 0.05;
dij.doseGrid.numOfVoxels = 1;
dij.ixDose{1} = 1;
w = 1;
scen = 1;
obj  = matRad_BEDProjection;
BED = computeSingleScenario(obj,dij,scen,w);
expectedResult = 2.4;

assertElementsAlmostEqual(BED,expectedResult,'absolute',1e-4);



function test_BED_computeSingleScenarioGradient
dij.physicalDose{1} = 2;
dij.mAlphaDose{1} = 1;
dij.mSqrtBetaDose{1} = sqrt(0.05)*2;
dij.ax{1} = 0.5;
dij.bx{1} = 0.05;
dij.doseGrid.numOfVoxels = 1;
dij.ixDose{1} = 1;
w = 1;
scen = 1;
obj  = matRad_BEDProjection;
doseGrad{1} = 1;
wGrad = projectSingleScenarioGradient(obj,dij,doseGrad,scen,w);

expectedResult = 2.8;

assertElementsAlmostEqual(wGrad,expectedResult,'absolute',1e-4);


function test_BED_setBiologicalDosePrescription
opti = DoseObjectives.matRad_SquaredOverdosing(1,2);
alphaX = 0.5;
betaX = 0.05;
opti = matRad_BEDProjection.setBiologicalDosePrescriptions(opti,alphaX,betaX);
assertEqual(opti.parameters{1},2.4 );

    
