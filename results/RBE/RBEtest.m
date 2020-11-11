clear

% Patient Data Import
load('BOXPHANTOM_LUNG_LARGE.mat');

% Treatment Plan

pln.radiationMode   = 'protons';     % either photons / protons / carbon
% pln.machine         = 'HITgenericRIFI3MMTOPAS_Water_forNiklas';
pln.machine         = 'Generic';

% modelName           = 'LEM';
modelName           = 'MCN';
quantityOpt         = 'RBExD';

pln.propOpt.bioOptimization = 'none';

pln.numOfFractions = 1;

pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 1;
% pln.propStf.longitudinalSpotSpacing = 10;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen'); % optimize on the nominal scenario

pln.heterogeneity.calcHetero = false;
pln.heterogeneity.useOriginalDepths = false;
pln.propDoseCalc.airOffsetCorrection = false;
%
stf = matRad_generateStfPencilBeam(pln,ct);
% stf = matRad_generateStf(ct,cst,pln);

%%
% dijMatRad = matRad_calcParticleDose(ct,stf,pln,cst);

%%
resultGUI_matRad = matRad_calcDoseDirect(ct,stf,pln,cst,1);

%% TOPAS heterogeneity correction
pln.propMC.proton_engine = 'TOPAS';
resultGUI_TOPAS = matRad_calcDoseDirectMC(ct,stf,pln,cst,0.001,1e2);
%%
% resultGUI_TOPAS_RBE = matRad_calcDoseDirectMC(ct,stf,pln,cst,weights,1e5);
%%
% dijTOPAS = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst);
% resultGUI_TOPAS = matRad_fluenceOptimization(dijTOPAS,cst,pln);