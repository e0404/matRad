function test_suite = test_generateBrachyStf

test_functions=localfunctions();

initTestSuite;

function pln = createPln(machine)
    % set up an example pln with relevant fields 
    pln.radiationMode   = 'brachy';
    pln.machine         = machine;
    
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

function test_generate_HDR()
    pln = createPln('HDR');
    stf = matRad_generateStf(ct, cst, pln, 0);
    assertTrue(isfield(stf, 'radiationMode'));
    assertTrue(isfield(stf, 'numOfSeedsPerNeedle'));
    assertTrue(isfield(stf, 'numOfNeedles'));
    assertTrue(isfield(stf, 'totalNumOfBixels'));
    assertTrue(isfield(stf, 'template'));
    assertTrue(isfield(stf, 'seedPoints'));

function test_generate_LDR()
    pln = createPln('LDR');
    stf = matRad_generateStf(ct, cst, pln, 0);
    assertTrue(isfield(stf, 'radiationMode'));
    assertTrue(isfield(stf, 'numOfSeedsPerNeedle'));
    assertTrue(isfield(stf, 'numOfNeedles'));
    assertTrue(isfield(stf, 'totalNumOfBixels'));
    assertTrue(isfield(stf, 'template'));
    assertTrue(isfield(stf, 'seedPoints'));


