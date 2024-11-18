clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

load('BOXPHANTOM.mat');

pln.radiationMode   = 'protons';
pln.machine         = 'Generic';
pln.multScen        = 'nomScen';
pln.numOfFractions  = 30;

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
pln.propDoseCalc.engine = 'HongPB';

%Optimization Settings
pln.propOpt.quantityOpt = 'RBExDose';

%% stf
stf = matRad_generateStf(ct,cst,pln);

%% Dose calc

% Select biological model
% Available models are:
%   RBEminMax (LET based): MCN, WED, CAR, LSM, (protons)
%                          HEL                 (helium)
%   kernel based:          LEM                 (carbon)

% We will compare the MCN model to constRBE.
% We will plan with the the MCN model.
pln.bioModel = matRad_MCNamara(); 
%alternative: pln.bioModel = matRad_bioModel(pln.radiationMode,'MCN');
%altnerative: pln.bioModel = 'MCN';
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% Fluence optimization
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Now let's recalculate with the constRBE model
pln.bioModel = matRad_ConstantRBE();
pln.bioModel.RBE = 1.1; %1.1 is standard, this is for illustration

resultGUI_recalc = matRad_calcDoseForward(ct,cst,stf,pln,resultGUI.w);

%% Compare Dose distributions
matRad_compareDose(resultGUI.RBExDose,resultGUI_recalc.RBExDose,ct,cst);