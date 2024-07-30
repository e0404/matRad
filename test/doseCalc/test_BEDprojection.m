function test_suite = test_BEDprojection

test_functions = localfunctions();

initTestSuite;

function test_BEDprojectionConstructor
bed = matRad_BEDProjection;
assertTrue(isobject(bed));
assertTrue(isa(bed, 'matRad_BEDProjection'));

%test compute single scenario BED projection 
function test_BED_computeSingleScenario
% ct of 3x3 with center voxel as target ,
% single pencil beam
a = sparse(zeros(9,1));
a(5) = 2;
dij.physicalDose{1} = a;
a(5)  = 1;
dij.mAlphaDose{1} = a;
a(5) = sqrt(0.05)*2;
dij.mSqrtBetaDose{1} = a;
dij.ax{1} = 0.5*ones(numel(a),1);
dij.bx{1} = 0.05*ones(numel(a),1);
dij.doseGrid.numOfVoxels = 9;
dij.ixDose{1} = 5;
w = 1;
scen = 1;
obj  = matRad_BEDProjection;
BED = computeSingleScenario(obj,dij,scen,w);
a = zeros(dij.doseGrid.numOfVoxels,1);
a(dij.ixDose{scen}) = 2.4;
expectedResult = a;

assertElementsAlmostEqual(BED,expectedResult,'absolute',1e-4);



function test_BED_computeSingleScenarioGradient
% ct of 3x3 with center voxel as target ,
% single pencil beam
a = sparse(zeros(9,1)); 
a(5) = 2;
dij.physicalDose{1} = a;
a(5)  = 1;
dij.mAlphaDose{1} = a;
a(5) = sqrt(0.05)*2;
dij.mSqrtBetaDose{1} = a;
dij.ax{1} = 0.5*ones(numel(a),1);
dij.bx{1} = 0.05*ones(numel(a),1);
dij.doseGrid.numOfVoxels = 9;
dij.ixDose{1} = 5;
w = 1;
scen = 1;
obj  = matRad_BEDProjection;
a = zeros(dij.doseGrid.numOfVoxels,1);
a(dij.ixDose{scen}) = 1;
doseGrad{1} = a ;
wGrad = projectSingleScenarioGradient(obj,dij,doseGrad,scen,w);

expectedResult = 2.8;

assertElementsAlmostEqual(wGrad,expectedResult,'absolute',1e-4);


function test_BED_setBiologicalDosePrescription
opti = DoseObjectives.matRad_SquaredOverdosing(1,2);
alphaX = 0.5;
betaX = 0.05;
opti = matRad_BEDProjection.setBiologicalDosePrescriptions(opti,alphaX,betaX);
assertEqual(opti.parameters{1},2.4 );

    
