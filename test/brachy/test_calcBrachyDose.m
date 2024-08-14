function test_suite = test_calcBrachyDose

test_functions =localfunctions();

initTestSuite;

% test if dij struct has required shape
function test_rightOutput()
    engine = DoseEngines.matRad_TG43BrachyEngine;
    pln.bioModel = matRad_bioModel('brachy', 'none');
    engine.assignPropertiesFromPln(pln);

    load PROSTATE.mat ct cst;
    load examplePlnAndStf.mat pln stf;
    pln.propDoseCalc.durationImplanted = Inf;
    
    
    dij = engine.setCalcDose(ct,cst,stf);
    assertTrue(isfield(dij,'doseGrid'));
    assertTrue(isfield(dij,'physicalDose'));
    assertTrue(iscell(dij.physicalDose));
    
    clear ct cst pln;
function test_getDistMtx
    % check calculated distance matrix against values calculated
    % on a spreadsheet (same folder as test) for 3x3 test point mtx
    seedPoints.x = [2,1,2];
    seedPoints.y = [2,3,1];
    seedPoints.z = [1,5,1];
    dosePoints.x = [0,1,1];
    dosePoints.y = [0,2,3];
    dosePoints.z = [0,1,2];
    verifyDistanceMatrix = [3,5.9,2.4;1,4.1,1.4;1.7,3,2.45];
    [funcDistanceMatrix,~] = matRad_getDistanceMatrix(seedPoints,dosePoints);
    
    assertElementsAlmostEqual(funcDistanceMatrix.dist,verifyDistanceMatrix,'absolute',0.1);

function test_getThetMtx
    engine = DoseEngines.matRad_TG43BrachyEngine;
    % check calculated angle matrix against values calculated
    % on a spreadsheet (same folder as test) for 3x3 test point mtx
    seedPoints.x = [2,1,2];
    seedPoints.y = [2,3,1];
    seedPoints.z = [1,5,1];
    dosePoints.x = [0,1,1];
    dosePoints.y = [0,2,3];
    dosePoints.z = [0,1,2];
    SeedDirection = sqrt(1/3)*ones(1,3);
    verifyThetaMatrix = [164.2,151.4,160.5;125.3,134.4,90;70.5,125.3,61.9];
    [DistanceMatrix,~] = matRad_getDistanceMatrix(seedPoints,dosePoints);
    [funcThetaMatrix,~] = engine.getThetaMatrix(SeedDirection,DistanceMatrix);
    
    assertElementsAlmostEqual(funcThetaMatrix,verifyThetaMatrix,'absolute',0.1)

function test_verifyDose1D

    engine = DoseEngines.matRad_TG43BrachyEngine;
    % verify if the calculated dose fits precalculated TG43.
    % The criterion is more than 99% gamma pass rate for 3mm and
    % 3%. For more information on gamma analysis, see help of the
    % matRad_gammaPassRate3D function
    
    % as a test, 10% of the points in a 81 x 81 grid are evaluated
    load referenceDoseCalculation refDos;
    
    machine = refDos.TG43_1D.basedata;
    r = refDos.coords.r;
    grids = refDos.coords.grids;
    
    doseCal = engine.getDoseRate1D_poly(machine,r);
    doseRef = refDos.TG43_1D.fullDose;
    
    [gammaPassRate,~] = ...
        matRad_gammaPassRate2D(doseRef,doseCal, grids, 3 , 3 , 0.1);
    
    assertTrue(gammaPassRate > 0.99)

function test_verifyDose2D

    engine = DoseEngines.matRad_TG43BrachyEngine;
    load referenceDoseCalculation refDos;
    
    machine = refDos.TG43_2D.basedata;
    grids = refDos.coords.grids;
    r = refDos.coords.r;
    theta = refDos.coords.theta;
    
    doseCal = engine.getDoseRate2D_poly(machine,r,theta);
    doseRef = refDos.TG43_2D.fullDose;
    
    [gammaPassRate,~] = ...
        matRad_gammaPassRate2D(doseRef,doseCal, grids, 3 , 3 , 0.1);
    
    assertTrue(gammaPassRate > 0.99)




