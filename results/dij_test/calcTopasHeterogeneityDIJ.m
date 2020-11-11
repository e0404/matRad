clear

% Patient Data Import
load('BOXPHANTOM_LUNG_LARGE_2e-1.mat');

% Treatment Plan

pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.machine         = 'generic_TOPAS_cropped_APM';
% pln.machine         = 'GenericAPM';

%
modelName           = 'none';
quantityOpt         = 'physicalDose';

pln.propOpt.bioOptimization = 'none';

%
pln.numOfFractions = 1;

pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 50;
pln.propStf.longitudinalSpotSpacing = 10;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen'); % optimize on the nominal scenario

pln.heterogeneity.calcHetero = true;
pln.heterogeneity.useOriginalDepths = false;
%
% stf = matRad_generateStfPencilBeam(pln,ct);
stf = matRad_generateStf(ct,cst,pln);

% Analytical heterogeneity correction
% weights = ones(1,sum([stf(:).totalNumOfBixels]));
% resultGUI_matRad = matRad_calcDoseDirect(ct,stf,pln,cst,weights);

%%
dij = matRad_calcParticleDose(ct,stf,pln,cst);
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%%
% cstHetero = matRad_cstHeteroAutoassign(cst);

% resultGUI_matRad_hetero = matRad_calcDoseDirect(ct,stf,pln,cstHetero,resultGUI_matRad.w);

%% TOPAS heterogeneity correction
% pln.propMC.proton_engine = 'TOPAS';
% resultGUI_TOPAS = matRad_calcDoseDirectMC(ct,stf,pln,cst,weights,1e5);
%%
dijTOPAS = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst);
resultGUI_TOPAS = matRad_fluenceOptimization(dijTOPAS,cst,pln);
%%
load('TOPASdijTestManyParticles.mat')
figure, plot(matRad_calcIDD(resultGUI.physicalDose)), hold on, plot(matRad_calcIDD(resultGUI_TOPAS.physicalDose))
legend({'matRad','TOPAS'},'Location','northwest')
ylabel('IDD dose')
save IDD.fig
figure, plot(matRad_calcIDD(resultGUI.physicalDose,'y',1)), hold on, plot(matRad_calcIDD(resultGUI_TOPAS.physicalDose,'y',1))
legend({'matRad','TOPAS'},'Location','northwest')
ylabel('profile dose')
save profiles.fig
figure, imagesc(resultGUI_TOPAS.physicalDose(:,:,80))
title('TOPAS dose')
save topasDose.fig
figure, imagesc(resultGUI.physicalDose(:,:,80))
title('matRad dose')
save matRadDose.fig
%%

% numOfSamples = 100;
% for i = 1:numOfSamples
%     ct_mod = matRad_modulateDensity(ct,cst,800);
%     i
%     resultGUI_mod{i} = matRad_calcDoseDirect(ct_mod,stf,pln,cst,resultGUI_matRad.w);                                                                      
% end
% 
% resultGUI_MC_mod.physicalDose = zeros(160,160,160);
% for i = 1:numOfSamples
%     resultGUI_MC_mod.physicalDose = resultGUI_MC_mod.physicalDose + resultGUI_mod{i}.physicalDose/numOfSamples;
% end
%%
% load CalcHomogeneousLungPhantom
% figure, plot(matRad_calcIDD(resultGUI_matRad.physicalDose)), hold on, plot(matRad_calcIDD(resultGUI_matRad_hetero.physicalDose))
% plot(matRad_calcIDD(resultGUI_TOPAS.physicalDose))
% plot(matRad_calcIDD(resultGUI_MC_mod.physicalDose))
% title('Homogeneous Phantom')
% legend({'matRad','matRad analytisch korrigiert','TOPAS','matRad modulated'},'Location','northwest')
% xlim([50 90])
% %%
% load CalcHeterogeneousLungPhantom
% figure, plot(matRad_calcIDD(resultGUI_matRad.physicalDose)), hold on, plot(matRad_calcIDD(resultGUI_matRad_hetero.physicalDose))
% plot(matRad_calcIDD(resultGUI_TOPAS.physicalDose))
% plot(matRad_calcIDD(resultGUI_MC_mod.physicalDose))
% title('Heterogeneous Phantom')
% legend({'matRad','matRad analytisch korrigiert','TOPAS','matRad modulated'},'Location','northwest')
% xlim([50 90])


