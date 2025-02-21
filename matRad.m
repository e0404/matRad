%% set matRad runtime configuration
matRad_rc

%% load patient data, e.g. TG119.mat
load BOXPHANTOM.mat

pln.numOfFractions  = 1;
pln.radiationMode   = 'VHEE';
pln.machine         = 'FermiEyges';
%pln.machine         = 'Focused';


pln.bioModel = 'none';
%pln.multScen = 'nomScen';

% Beam geometry settings
%pln.propStf.generator = 'ParticleVHEE';
pln.propStf.energy  = 200;  % desired energy in MeV
pln.propStf.bixelWidth   = 5;
pln.propStf.gantryAngles = [0];
pln.propStf.gantryAngles = [0:72:359];
pln.propStf.couchAngles  = [0,0,0,0,0];
pln.propStf.numOfBeams   = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter    = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

% Dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5;
pln.propDoseCalc.doseGrid.resolution.y = 5;
pln.propDoseCalc.doseGrid.resolution.z = 5;
%% 

% Optimization settings
pln.propOpt.quantityOpt  = 'physicalDose';
pln.propOpt.optimizer    = 'IPOPT';
pln.propOpt.runDAO       = false;
pln.propSeq.runSequencing = false;

%% Generate steering file and calculate dose
stf = matRad_generateStf(ct,cst,pln);
%%
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
%%
resultGUI = matRad_fluenceOptimization(dij, cst, pln);
%%
%resultGUI = matRad_sequencing(resultGUI, stf, dij, pln);

% %% (Optionally) Run DAO and GUI
% if strcmp(pln.radiationMode, 'photons') && pln.propOpt.runDAO
%    resultGUI = matRad_directApertureOptimization(dij, cst, resultGUI.apertureInfo, resultGUI, pln);
%    matRad_visApertureInfo(resultGUI.apertureInfo);
% end

%matRadGUI;
% resultGUI = matRad_planAnalysis(resultGUI, ct, cst, stf, pln);