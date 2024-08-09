matRad_rc;
matRad_cfg = MatRad_Config.instance();

load('BOXPHANTOM.mat');

pln.radiationMode = 'protons';
pln.machine       = 'Generic';
pln.propDoseCalc.calcLET = 0;

pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);

pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 8;
pln.propDoseCalc.doseGrid.resolution.y = 8;
pln.propDoseCalc.doseGrid.resolution.z = 8;

pln.multScen = matRad_multScen(ct, 'nomScen');

%% stf
stf = matRad_generateStf(ct,cst,pln);

%% Dose calc

% Select biological model
% Available models are:
%   RBEminMax (LET based): MCN, WED, CAR, LSM, (protons)
%                          HEL                 (helium)
%   kernel based:          LEM                 (carbon)
%   Tabulated              TAB                  (proton, helium, carbon) (requires RBEtable and spectra in base data)

% Example 1: load the MCN model. The machine input is optional in this case.
% The output is an instance of the model. 
pln.bioModel = matRad_bioModel(pln.radiationMode,'MCN', pln.machine);

% Alternatively, the biological model can also be assigned as a structure.
% The only andatory field in this case is the 'model' field:
% pln.bioModel.model = 'MCN';

% Example 2: Some models have parameters that can be tuned by the user.
% For example, we can instantiate a constRBE model
%
% pln.bioModel = matRad_bioModel(pln.radiationMode, 'constRBE');
%
% and assign a custom value for the constant RBE
%
% pln.bioModel.RBE = 4.2;


% Example 3: Tabulated models require  the additional specification of an
% RBE table and additional base data entries. Additional details are
% provided in the Wiki

% pln.bioModel = matRad_bioModel(pln.radiationMode, 'TAB', pln.machine);

% pln.bioModel.RBEtableName       = 'RBEtable_rapidLEMI_testTable';
% pln.bioModel.fragmentsToInclude  = {'H'};


pln.propDoseCalc.engine = 'HongPB';
dij_PB = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% Fluence opt
pln.propOpt.quantityOpt = 'RBExD';
resultGUI = matRad_fluenceOptimization(dij_PB,cst,pln);


matRadGUI;