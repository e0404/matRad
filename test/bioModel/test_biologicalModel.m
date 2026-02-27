function test_suite = test_biologicalModel

test_functions = localfunctions();

initTestSuite;

% function test_bioModelConstructor
%     if moxunit_util_platform_is_octave()
%          assertExceptionThrown(@() matRad_BiologicalModel());
%          assertExceptionThrown(@() matRad_BiologicalModel('photons'));
%     else
%         assertExceptionThrown(@() matRad_BiologicalModel(),'MATLAB:minrhs');
%         assertExceptionThrown(@() matRad_BiologicalModel('photons'),'MATLAB:minrhs');
%     end


function test_setBiologicalModel
    bioModel = matRad_bioModel('photons','none');
    assertTrue(isa(bioModel, 'matRad_EmptyBiologicalModel'));
    bioModel = matRad_bioModel('protons', 'none');
    assertTrue(isa(bioModel, 'matRad_EmptyBiologicalModel'));
    
    if moxunit_util_platform_is_octave()
        assertExceptionThrown(@(model) matRad_bioModel('photons', 'MCN'));
        assertExceptionThrown(@(model) matRad_bioModel('protons', 'HEL'));
    else
        assertExceptionThrown(@(model) matRad_bioModel('photons', 'MCN'), 'matRad:Error');
        assertExceptionThrown(@(model) matRad_bioModel('protons', 'HEL'),'matRad:Error');
    end

function test_setBiologicalModelProvidedQuantities
    bioModel = matRad_bioModel('protons', 'MCN', {'physicalDose','LET'});
    assertTrue(isa(bioModel, 'matRad_MCNamara'));
    
    if moxunit_util_platform_is_octave()
        assertExceptionThrown(@(model) matRad_bioModel('protons', 'MCN', {'physicalDose'}));        
    else
        assertExceptionThrown(@(model) matRad_bioModel('photons', 'MCN', {'physicalDose'}), 'matRad:Error');
    end

function test_tissueParameters_emptyModel
    bioModel = matRad_EmptyBiologicalModel();
    abx = bioModel.getAvailableTissueParameters(struct());
    assertTrue(isempty(abx));

function test_tissueParameters_kernelModel
    bioModel = matRad_KernelBasedLEM();
    abx = bioModel.getAvailableTissueParameters(struct('machine','Generic','radiationMode','carbon'));
    assertTrue(isnumeric(abx));
    assertEqual(size(abx,2),2);
    assertTrue(size(abx,1) >= 1);
    
    if moxunit_util_platform_is_octave()
         assertExceptionThrown(@() bioModel.getAvailableTissueParameters(struct('machine','Generic','radiationMode','photons')));
    else
        assertExceptionThrown(@() bioModel.getAvailableTissueParameters(struct('machine','Generic','radiationMode','photons')), 'matRad:Error');
    end


function test_calcBiologicalQuantitiesForBixel_MCN
    bioModel = matRad_bioModel('protons','MCN');
    kernels.LET = 0.9810;
        
    bixel.energyIx = 11;
    load protons_Generic.mat;
    bixel.baseData = machine.data(11);
    bixel.radDepths = 0;
    bixel.vAlphaX = 0.5;
    bixel.vBetaX = 0.05;
    bixel.vTissueIndex = 1;
    
    test_bixelAlpha = 0.5170;
    test_bixelBeta = 0.0593;
    
    [bixel] = calcBiologicalQuantitiesForBixel(bioModel,bixel,kernels);
    assertTrue(isstruct(bixel));
    assertTrue(isnumeric(bixel.alpha));
    assertTrue(isnumeric(bixel.beta));
    assertElementsAlmostEqual(bixel.alpha,test_bixelAlpha,'absolute', 1e-4);
    assertElementsAlmostEqual(bixel.beta, test_bixelBeta,'absolute', 1e-4);

function test_calcBiologicalQuantitiesForBixel_WED
    bioModel = matRad_bioModel('protons','WED');
    kernels.LET = 0.9810;
        
    bixel.energyIx = 11;
    load protons_Generic.mat;
    bixel.baseData = machine.data(11);
    bixel.radDepths = 0;
    bixel.vAlphaX = 0.5;
    bixel.vBetaX = 0.05;
    bixel.vTissueIndex = 1;
    
    test_bixelAlpha = 0.5213;
    test_bixelBeta = 0.0500;
    
    [bixel] = calcBiologicalQuantitiesForBixel(bioModel,bixel,kernels);
    assertTrue(isstruct(bixel));
    assertTrue(isnumeric(bixel.alpha));
    assertTrue(isnumeric(bixel.beta));
    assertElementsAlmostEqual(bixel.alpha,test_bixelAlpha,'absolute', 1e-4);
    assertElementsAlmostEqual(bixel.beta, test_bixelBeta,'absolute', 1e-4);

function test_calcBiologicalQuantitiesForBixel_HEL
    bioModel = matRad_bioModel('helium','HEL');
    kernels.LET = 4.1914;
        
    bixel.energyIx = 11;
    load helium_Generic.mat;
    bixel.baseData = machine.data(11);
    bixel.radDepths = 0;
    bixel.vAlphaX = 0.5;
    bixel.vBetaX = 0.05;
    bixel.vTissueIndex = 1;
    
    test_bixelAlpha = 0.5190;
    test_bixelBeta = 0.0500;
    
    [bixel] = calcBiologicalQuantitiesForBixel(bioModel,bixel,kernels);
    assertTrue(isstruct(bixel));
    assertTrue(isnumeric(bixel.alpha));
    assertTrue(isnumeric(bixel.beta));
    assertElementsAlmostEqual(bixel.alpha,test_bixelAlpha,'absolute', 1e-4);
    assertElementsAlmostEqual(bixel.beta, test_bixelBeta,'absolute', 1e-4);

function test_calcBiologicalQuantitiesForBixel_LEM
    bioModel = matRad_bioModel('carbon','LEM');
    kernels.LET = 34.968;
    kernels.alpha = [0.2978 0.7053];
    kernels.beta = [0.0412 0.0443];

    bixel.energyIx = 4;
    load carbon_Generic.mat;
    bixel.baseData = machine.data(4);
    bixel.radDepths = 0;
    bixel.vAlphaX = 0.5;
    bixel.vBetaX = 0.05;
    bixel.vTissueIndex = 1;
    
    test_bixelAlpha = 0.2978;
    test_bixelBeta = 0.0412;
    
    [bixel] = calcBiologicalQuantitiesForBixel(bioModel,bixel,kernels);
    assertTrue(isstruct(bixel));
    assertTrue(isnumeric(bixel.alpha));
    assertTrue(isnumeric(bixel.beta));
    assertElementsAlmostEqual(bixel.alpha,test_bixelAlpha,'absolute', 1e-4);
    assertElementsAlmostEqual(bixel.beta, test_bixelBeta,'absolute', 1e-4);

% function test_bioOptimization_MCN_BED
%     matRad_rc;
%     matRad_cfg = MatRad_Config.instance();
% 
%     load('BOXPHANTOM.mat');
% 
%     pln.radiationMode = 'protons';
%     pln.machine       = 'Generic';
% 
%     pln.numOfFractions        = 30;
%     pln.propStf.gantryAngles  = 0;
%     pln.propStf.couchAngles   = 0;
%     pln.propStf.bixelWidth    = 5;
%     pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
% 
%     pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
%     pln.propOpt.runDAO        = 0;
%     pln.propSeq.runSequencing = 0;
% 
%     % dose calculation settings
%     pln.propDoseCalc.doseGrid.resolution.x = 8;
%     pln.propDoseCalc.doseGrid.resolution.y = 8;
%     pln.propDoseCalc.doseGrid.resolution.z = 8;
%     pln.propDoseCalc.calcLET = false;
% 
%     pln.multScen = matRad_multScen(ct, 'nomScen');
% 
%     stf = matRad_generateStf(ct,cst,pln);
%     %Dose calc
%     pln.bioModel = matRad_bioModel(pln.radiationMode,'MCN');
% 
%     pln.propDoseCalc.engine = 'HongPB';
%     dij = matRad_calcDoseInfluence(ct,cst,stf,pln);
% 
%     %Opt
%     matRad_cfg.defaults.propOpt.maxIter = 10;
% 
%     pln.propOpt.quantityOpt = 'BED';
%     resultGUI = matRad_fluenceOptimization(dij,cst,pln);
%     cumulativeBED = sum(resultGUI.BED, 'all');
% 
%     expectedCumulativeBED = 1.271963512370974e+05;
% 
%     assertTrue(isstruct(resultGUI));
%     assertTrue(isfield(resultGUI, 'BED'));
%     assertElementsAlmostEqual(cumulativeBED,expectedCumulativeBED, 1e-04);