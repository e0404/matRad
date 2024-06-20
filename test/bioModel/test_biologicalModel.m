function test_suite = test_biologicalModel

test_functions = localfunctions();

initTestSuite;

function test_bioModelConstructor
    if moxunit_util_platform_is_octave()
         assertExceptionThrown(@() matRad_BiologicalModel());
         assertExceptionThrown(@() matRad_BiologicalModel('photons'));
    else
        assertExceptionThrown(@() matRad_BiologicalModel(),'MATLAB:minrhs');
        assertExceptionThrown(@() matRad_BiologicalModel('photons'),'MATLAB:minrhs');
    end


function test_setBiologicalModel 
    bioModel = matRad_BiologicalModel('photons','RBExD','LEM');
    assertTrue(strcmp(bioModel.model, 'none')); %default photon model
    bioModel = matRad_BiologicalModel('protons','RBExD','LEM');
    assertTrue(strcmp(bioModel.model, 'constRBE')); %default proton model
    bioModel = matRad_BiologicalModel('carbon','RBExD','MCN');
    assertTrue(strcmp(bioModel.model, 'LEM')); %default carbon model

function test_calcLQParameter
    bioModel = matRad_BiologicalModel('photons','RBExD','LEM');
    if moxunit_util_platform_is_octave()
        assertExceptionThrown(@() bioModel.calcLQParameter());
        assertExceptionThrown(@() bioModel.calcLQParameterForKernel());
    else
        assertExceptionThrown(@() bioModel.calcLQParameter(),'MATLAB:minrhs');
        assertExceptionThrown(@() bioModel.calcLQParameterForKernel(),'MATLAB:minrhs');
    end
    

function test_calcLQParameterForKernel_protonsMCN
    bioModel = matRad_BiologicalModel('protons','RBExD','MCN');
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
    
    [bixelAlpha,bixelBeta] = calcLQParameterForKernel(bioModel,bixel,kernels);
    assertTrue(isnumeric(bixelBeta));
    assertTrue(isnumeric(bixelAlpha));
    assertElementsAlmostEqual(bixelAlpha,test_bixelAlpha,'absolute', 1e-4);
    assertElementsAlmostEqual(bixelBeta, test_bixelBeta,'absolute', 1e-4);


function test_calcLQParameterForKernel_protonsWED
    bioModel = matRad_BiologicalModel('protons','RBExD','WED');
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
    
    [bixelAlpha,bixelBeta] = calcLQParameterForKernel(bioModel,bixel,kernels);
    
    assertTrue(isnumeric(bixelBeta));
    assertTrue(isnumeric(bixelAlpha));
    assertElementsAlmostEqual(bixelAlpha,test_bixelAlpha,'absolute', 1e-4);
    assertElementsAlmostEqual(bixelBeta, test_bixelBeta,'absolute', 1e-4);

function test_calcLQParameterForKernel_heliumHEL
    bioModel = matRad_BiologicalModel('helium','RBExD','HEL');
    kernels.LET = 4.1914;
        
    bixel.energyIx = 11;
    load helium_Generic.mat;
    bixel.baseData = machine.data(11);
    bixel.radDepths = 0.05;
    bixel.vAlphaX = 0.5;
    bixel.vBetaX = 0.05;
    bixel.vTissueIndex = 1;
    
    test_bixelAlpha = 0.5190;
    test_bixelBeta = 0.0500;
    
    [bixelAlpha,bixelBeta] = calcLQParameterForKernel(bioModel,bixel,kernels);
    assertTrue(isnumeric(bixelBeta));
    assertTrue(isnumeric(bixelAlpha));
    assertElementsAlmostEqual(bixelAlpha,test_bixelAlpha,'absolute', 1e-4);
    assertElementsAlmostEqual(bixelBeta, test_bixelBeta,'absolute', 1e-4);

function test_calcLQParameterForKernel_carbonLEM
    bioModel = matRad_BiologicalModel('carbon','RBExD','LEM');
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
    
    [bixelAlpha,bixelBeta] = calcLQParameterForKernel(bioModel,bixel,kernels);
    assertTrue(isnumeric(bixelBeta));
    assertTrue(isnumeric(bixelAlpha));
    assertElementsAlmostEqual(bixelAlpha,test_bixelAlpha,'absolute', 1e-4);
    assertElementsAlmostEqual(bixelBeta, test_bixelBeta,'absolute', 1e-4);

