function test_suite = test_Analytical

test_functions=localfunctions();

initTestSuite;

function test_getAnalyticalEngineFromPln
    % Single gaussian lateral model
    protonDummyPln = struct('radiationMode','protons','machine','Generic');
    protonDummyPln.propDoseCalc.engine = 'AnalyticalPB';
    engine = DoseEngines.matRad_ParticleAnalyticalBortfeldEngine.getEngineFromPln(protonDummyPln);
    assertTrue(isa(engine,'DoseEngines.matRad_ParticleAnalyticalBortfeldEngine'));

    % Double Gaussian lateral model 
    % If you don't have my clusterDose basedata you cannot try this :P
    %{
    protonDummyPln = struct('radiationMode','protons','machine','Generic_clusterDose');
    protonDummyPln.propDoseCalc.engine = 'AnalyticalPB';
    engine = DoseEngines.matRad_ParticleAnalyticalBortfeldEngine.getEngineFromPln(protonDummyPln);
    assertTrue(isa(engine,'DoseEngines.matRad_ParticleAnalyticalBortfeldEngine'));
    %}

function test_loadMachineForAnalytical
    possibleRadModes = DoseEngines.matRad_ParticleAnalyticalBortfeldEngine.possibleRadiationModes;
    for i = 1:numel(possibleRadModes)
        machine = DoseEngines.matRad_ParticleAnalyticalBortfeldEngine.loadMachine(possibleRadModes{i},'Generic');
        assertTrue(isstruct(machine));
        assertTrue(isfield(machine, 'meta'));
        assertTrue(isfield(machine.meta, 'radiationMode'));
        assertTrue(strcmp(machine.meta.radiationMode, 'protons'));
    end

function test_calcDoseAnalytical
    matRad_cfg = MatRad_Config.instance();

    protonDummyPln = struct('radiationMode','protons','machine','Generic');
    protonDummyPln.propDoseCalc.engine = 'AnalyticalPB';

    load([protonDummyPln.radiationMode '_' protonDummyPln.machine]);

    load BOXPHANTOM.mat

    stf = matRad_generateStf(ct, cst, protonDummyPln);
    
    resultGUI = matRad_calcDoseForward(ct, cst, stf, protonDummyPln, ones([1, stf(:).totalNumOfBixels]));

    assertTrue(isfield(resultGUI, 'physicalDose'));
    assertTrue(isfield(resultGUI, 'w'));
    assertTrue(isequal(size(ct.cube{1}), size(resultGUI.physicalDose)))

function test_nonSupportedSettings
    % Radiation mode other than protons not implemented 
    carbonDummyPln = struct('radiationMode','carbon','machine','Generic');
    carbonDummyPln.propDoseCalc.engine = 'AnalyticalPB';
    engine = DoseEngines.matRad_ParticleAnalyticalBortfeldEngine.getEngineFromPln(carbonDummyPln);
    assertTrue(~isa(engine,'DoseEngines.matRad_ParticleAnalyticalBortfeldEngine'));

    % Biological models, LET, other lateral models not implemented
    protonDummyPln = struct('radiationMode','protons','machine','Generic');
    protonDummyPln.propDoseCalc.engine = 'AnalyticalPB';
    protonDummyPln.propDoseCalc.calcLET = true;
    protonDummyPln.propDoseCalc.calcBioDose = true;
    protonDummyPln.propDoseCalc.lateralModel = 'double';
    engine = DoseEngines.matRad_ParticleAnalyticalBortfeldEngine.getEngineFromPln(protonDummyPln);
    assertTrue(isa(engine,'DoseEngines.matRad_ParticleAnalyticalBortfeldEngine'));
    load BOXPHANTOM.mat
    stf = matRad_generateStf(ct, cst, protonDummyPln);
    resultGUI = matRad_calcDoseForward(ct, cst, stf, protonDummyPln, ones([1, stf(:).totalNumOfBixels]));
    assertTrue(~engine.calcLET)
    %assertTrue(~engine.calcBioDose) % Access protected property

    % Invalid machine without radiation mode field
    matRad_cfg = MatRad_Config.instance();
    protonDummyPln = struct('radiationMode','protons','machine','Empty');
    protonDummyPln.propDoseCalc.engine = 'AnalyticalPB';
    machine = [];
    assertExceptionThrown(@() DoseEngines.matRad_ParticleAnalyticalBortfeldEngine.isAvailable(protonDummyPln, machine));



