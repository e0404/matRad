% clear

% Patient Data Import
load('BOXPHANTOM_LUNG_LARGE.mat');

% Treatment Plan

pln.radiationMode   = 'protons';     % either photons / protons / carbon
% pln.machine         = 'HITgenericRIFI3MMTOPAS_Water_forNiklas';
pln.machine         = 'Generic';

% modelName           = 'LEM';
modelName           = 'MCN';
% quantityOpt         = 'RBExD';
quantityOpt         = 'physicalDose';


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
pln.propDoseCalc.airOffsetCorrection = true;
%
stf = matRad_generateStfPencilBeam(pln,ct);
% stf = matRad_generateStf(ct,cst,pln);
weights = ones(1,pln.propStf.numOfBeams);
% %%
% dijMatRad = matRad_calcParticleDose(ct,stf,pln,cst,1);
% dijTOPAS = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst,1e4,1);
%%
resultGUI_matRad = matRad_calcDoseDirect(ct,stf,pln,cst,weights);

%% TOPAS heterogeneity correction
pln.propMC.proton_engine = 'MCsquare';
resultGUI_MCsquare = matRad_calcDoseDirectMC(ct,stf,pln,cst,1,1e4);
%%
pln.propMC.proton_engine = 'TOPAS';
tic
resultGUI_TOPAS = matRad_calcDoseDirectMC(ct,stf,pln,cst,1,1e4);
toc
%%
% resultGUI_TOPAS_RBE = matRad_calcDoseDirectMC(ct,stf,pln,cst,weights,1e5);
%%
% dijTOPAS = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst);
% resultGUI_TOPAS = matRad_fluenceOptimization(dijTOPAS,cst,pln);

figure, plot(matRad_calcIDD(resultGUI_matRad.physicalDose)), hold on,...
    plot(matRad_calcIDD(resultGUI_MCsquare.physicalDose)),...
    plot(matRad_calcIDD(resultGUI_TOPAS.physicalDose))
legend({'matRad','MCsquare','TOPAS'},'Location','northwest')
xlim([10 100])
ylabel('proton physicalDose')
