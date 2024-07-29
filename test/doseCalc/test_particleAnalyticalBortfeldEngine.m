function test_suite = test_baseEngine

test_functions=localfunctions();

initTestSuite;

function test_properties
    %% Set matRad and pln
    matRad_rc

    %load HEAD_AND_NECK
    load TG119.mat
    %load PROSTATE.mat
    %load LIVER.mat
    %load BOXPHANTOM.mat

    % meta information for treatment plan
    pln.numOfFractions  = 30;
    pln.radiationMode   = 'protons';           % either photons / protons / helium / carbon
    pln.machine         = 'Generic';
    pln.propDoseCalc.engine = 'AnalyticalPB';

    % beam geometry settings
    pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
    pln.propStf.gantryAngles    = [0:72:359]; % [°] ;
    pln.propStf.couchAngles     = [0 0 0 0 0]; % [°] ;
    pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
    pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
    % optimization settings
    pln.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
    pln.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

    quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
    modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons
    % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons
    % LEM: Local Effect Model for carbon ions       % HEL: data-driven RBE parametrization for helium
    % dose calculation settings
    pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
    pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
    pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

    scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'

    % retrieve bio model parameters
    %pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

    % retrieve scenarios for dose calculation and optimziation
    pln.multScen = matRad_multScen(ct,scenGenType);

    % optimization settings
    pln.propOpt.optimizer       = 'IPOPT';
    % pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
    %                                       % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
    pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
    pln.propSeq.runSequencing   = true;  % true: run sequencing, false: don't / will be ignored for particles and also triggered by runDAO below

    %% Generate stf
    stf = matRad_generateStf(ct,cst,pln);

    %% Set LET true. Try to run a dose calculation to throw a warning
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
    engine.calcLET = true;

    %call the calcDose funktion
    %dij = engine.calcDoseInfluence(ct,cst,stf);

    %% Verify
    assertExceptionThrown(@() engine.calcDoseInfluence(ct,cst,stf), 'Engine does not support LET calculation! Disabling!');

%{
function test_doseEngineBaseAbstract
    if moxunit_util_platform_is_octave()
        assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase(),'');
    else
        assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase(),'MATLAB:class:abstract');
    end

function test_abstractIsAvailable
    assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase.isAvailable,'matRad:Error');

function test_getAvailableEngines
    avail = DoseEngines.matRad_DoseEngineBase.getAvailableEngines();
    assertTrue(~isempty(avail));
    assertTrue(isstruct(avail));
    assertTrue(iscolumn(avail));

    photonDummyPln = struct('radiationMode','photons','machine','Generic');
    availPhotons = DoseEngines.matRad_DoseEngineBase.getAvailableEngines(photonDummyPln);
    assertTrue(~isempty(avail));
    assertTrue(isstruct(avail));
    assertTrue(iscolumn(avail));
    assertTrue(numel(availPhotons) < numel(avail));

function test_loadMachine
    machine = DoseEngines.matRad_DoseEngineBase.loadMachine('photons','Generic');
    assertTrue(isstruct(machine));

    if moxunit_util_platform_is_octave()
        assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase.loadMachine());
        assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase.loadMachine('photons'));
    else
        assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase.loadMachine(),'MATLAB:minrhs');
        assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase.loadMachine('photons'),'MATLAB:minrhs');
    end
    assertExceptionThrown(@() DoseEngines.matRad_DoseEngineBase.loadMachine('grbl','grbl'),'matRad:Error');

function test_getEngineFromPlnDefaults
    photonDummyPln = struct('radiationMode','photons','machine','Generic');
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(photonDummyPln);
    assertTrue(isa(engine,'DoseEngines.matRad_PhotonPencilBeamSVDEngine'));

    protonDummyPln = struct('radiationMode','protons','machine','Generic');
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(protonDummyPln);
    assertTrue(isa(engine,'DoseEngines.matRad_ParticleHongPencilBeamEngine'));

    carbonDummyPln = struct('radiationMode','carbon','machine','Generic');
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(carbonDummyPln);
    assertTrue(isa(engine,'DoseEngines.matRad_ParticleHongPencilBeamEngine'));

    heliumDummyPln = struct('radiationMode','helium','machine','Generic');
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(heliumDummyPln);
    assertTrue(isa(engine,'DoseEngines.matRad_ParticleHongPencilBeamEngine'));

function test_getEngineFromPlnByName
    protonDummyPln = struct('radiationMode','protons','machine','Generic','propDoseCalc',struct('engine','MCsquare'));
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(protonDummyPln);
    assertTrue(isa(engine,'DoseEngines.matRad_ParticleMCsquareEngine'));

    %Wrong name
    photonDummyPln = struct('radiationMode','photons','machine','Generic','propDoseCalc',struct('engine','MCsquare'));
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(photonDummyPln);
    assertTrue(isa(engine,'DoseEngines.matRad_PhotonPencilBeamSVDEngine'));

%}


