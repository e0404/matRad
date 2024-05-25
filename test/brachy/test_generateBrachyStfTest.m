function test_suite = test_generateBrachyStfTest

test_functions=localfunctions();

initTestSuite;



function test_addBrachyInfo()
    
    global plnTest;
    

    % set up an example pln with relevant fields 
    pln.radiationMode   = 'brachy';
    pln.machine         = 'HDR';
    
    % geometry settings
    load PROSTATE.mat ct cst;
    pln.propStf.templateRoot             = matRad_getTemplateRoot(ct, cst);
    pln.propStf.needle.seedDistance      = 1; % [mm] seed distance on needle
    pln.propStf.needle.seedsNo           = 2; % number of seeds per needle
    pln.propStf.template.numOfXPoints    = 2;
    pln.propStf.template.numOfYPoints    = 2;
    pln.propStf.template.xScale          = 1; % [mm] distance of neighbouring points
    pln.propStf.template.yScale          = 1; % [mm] distance of neighbouring points
    pln.propStf.template.activeNeedles   = [0 0 0 0 0 0 0 0 0 0 0 0 0;...
                                            0 0 0 0 0 0 0 0 0 0 0 0 0;...
                                            0 0 0 0 1 1 0 1 1 0 0 0 0;...
                                            1 0 1 0 1 0 0 0 1 0 1 0 1;...
                                            0 1 0 1 0 1 0 1 0 1 0 1 0;...
                                            1 0 1 0 1 0 0 0 1 0 1 0 1;...
                                            0 1 0 1 0 1 0 1 0 1 0 1 0;...
                                            1 0 1 0 1 0 1 0 1 0 1 0 1;...
                                            0 1 0 1 0 1 0 1 0 1 0 1 0;...
                                            0 0 1 0 1 0 0 0 1 0 1 0 0;...
                                            0 0 0 1 0 0 0 0 0 1 0 0 0;...
                                            0 0 0 0 0 0 0 0 0 0 0 0 0;...                                                                      
                                            0 0 0 0 0 0 0 0 0 0 0 0 0];
    pln.propStf.bixelWidth = 10;
    plnTest = pln;



function test_rightOutput()
    global plnTest;
    load PROSTATE.mat ct cst;
    stf = matRad_generateStf(ct, cst, plnTest, 0);
    assertTrue(isfield(stf, 'radiationMode'));
    assertTrue(isfield(stf, 'numOfSeedsPerNeedle'));
    assertTrue(isfield(stf, 'numOfNeedles'));
    assertTrue(isfield(stf, 'totalNumOfBixels'));
    assertTrue(isfield(stf, 'template'));
    assertTrue(isfield(stf, 'seedPoints'));