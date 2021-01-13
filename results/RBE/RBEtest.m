clear

% Patient Data Import
load('BOXPHANTOM.mat');
% load('BOXPHANTOM_LUNG_LARGE')
% ct.cube{1} = ones(160,160,160);
% ct.cubeHU{1} = zeros(160,160,160);
% cst{1,4}{1,1} = [1:prod([160,160,160])]';
% Treatment Plan

pln.radiationMode   = 'carbon';
% pln.radiationMode   = 'helium';     % either photons / protons / carbon
% pln.machine         = 'HITfixedBL';
pln.machine         = 'HITgenericRIFI3MMTOPAS_Water_forNiklas';

% modelName           = 'none';
modelName           = 'LEM';
quantityOpt         = 'RBExD';
% quantityOpt         = 'physicalDose';


pln.propOpt.bioOptimization = 'none';

pln.numOfFractions = 15;

pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 7;
pln.propStf.longitudinalSpotSpacing = 1;
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
%  stf = matRad_generateStf(ct,cst,pln);
% 
%%
dijMatRad = matRad_calcParticleDose(ct,stf,pln,cst);
resultGUI_matRad = matRad_fluenceOptimization(dijMatRad,cst,pln);
% dijTOPAS = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst,1e4,1);
%%

% stf.ray.energy = machine.data(60).energy;
% weights = ones(1,stf.totalNumOfBixels);
% resultGUI_matRad_1 = matRad_calcDoseDirect(ct,stf,pln,cst,weights);

%% TOPAS heterogeneity correction
% if strcmp(pln.radiationMode,'protons')
%     pln.propMC.proton_engine = 'MCsquare';
%     resultGUI_MCsquare = matRad_calcDoseDirectMC(ct,stf,pln,cst,1,1e5);
% end
%%
% pln.propMC.carbon_engine = 'TOPAS';
tic
resultGUI_TOPAS_LEM1 = matRad_calcDoseDirectMC(ct,stf,pln,cst,resultGUI_matRad.w,1e6);
toc
%%
% save CarbonSOBP_LEM
%%
% figure, plot(matRad_calcIDD(resultGUI_TOPAS_LEM1_noSpec.physicalDose,'y',1),'--'), hold on, plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.physicalDose,'y',1),'--')

% figure, plot(matRad_calcIDD(resultGUI_TOPAS_LEM1_noSpec.RBExD,'y',1),'--'), hold on, plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.RBExD,'y',1),'--')

%% Profiles

figure, plot(matRad_calcIDD(resultGUI_matRad.physicalDose,'y',1),'b'), hold on, plot(matRad_calcIDD(resultGUI_matRad.RBExD,'y',1),'b--')
plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.physicalDose,'y',1),'m')
plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.RBExD,'y',1),'m--')
legend({'matRad physDose','matRad RBExD','TOPAS physDose','TOPAS RBExD'})

%% IDD
figure, plot(matRad_calcIDD(resultGUI_matRad.physicalDose),'b'), hold on, plot(matRad_calcIDD(resultGUI_matRad.RBExD),'b--')
plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.physicalDose),'m')
plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.RBExD),'m--')
legend({'matRad physDose','matRad RBExD','TOPAS physDose','TOPAS RBExD'})


%% alpha beta
figure, plot(matRad_calcIDD(resultGUI_matRad.alpha,'y',1),'b'), hold on, plot(sqrt(matRad_calcIDD(resultGUI_matRad.beta,'y',1)),'b--')
plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.alpha,'y',1),'m')
plot(sqrt(matRad_calcIDD(resultGUI_TOPAS_LEM1.beta,'y',1)),'m--')
legend({'matRad alpha','matRad sqrt(beta)','TOPAS alpha','TOPAS sqrt(beta)'})

%% RBE
figure, plot(matRad_calcIDD(resultGUI_matRad.RBE,'y',1),'b'), hold on, plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.RBE,'y',1),'m')
legend({'matRad RBE','TOPAS RBE'})


%% RBExD manuell berechnet
% 
% dij.ax = full(reshape(dijMatRad.ax,[160 160 160]));
% dij.bx = full(reshape(dijMatRad.bx,[160 160 160]));
% ix = dij.bx~=0;
% RBExD = zeros([160 160 160]);
% RBExD_TOPAS = zeros([160 160 160]);
% 
% RBExD(ix) =         (sqrt(dij.ax(ix).^2 + 4*dij.bx(ix) .* resultGUI_matRad.effect(ix))      -dij.ax(ix))./(2.*dij.bx(ix));
% RBExD_TOPAS(ix) =   (sqrt(dij.ax(ix).^2 + 4*dij.bx(ix) .* resultGUI_TOPAS_LEM1.effect(ix))  -dij.ax(ix))./(2.*dij.bx(ix));
% 
% figure, plot(matRad_calcIDD(resultGUI_matRad.RBExD,'y',1),'r'), hold on, plot(matRad_calcIDD(RBExD,'y',1),'c--'), plot(matRad_calcIDD(RBExD_TOPAS,'y',1),'b--')
% legend({'matRad RBExD','matRad RBExD manuell','TOPAS RBExD','TOPAS RBExD manuell'})
% figure%, hold on, plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.RBExD,'y',1))
% plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.physicalDose,'y',1),'--')
% plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.RBExD,'y',1))
% % plot(matRad_calcIDD(resultGUI_TOPAS_direct_old.physicalDose),'--')
% % plot(matRad_calcIDD(resultGUI_TOPAS_direct_old.RBExD))
% legend({'matRad physDose','matRad','TOPAS LEM1 physDose','TOPAS LEM1'})
% % legend({'matRad physDose','matRad','TOPAS amtrack physDose','TOPAS amtrack'})
% %%
% figure, plot(matRad_calcIDD(resultGUI_matRad.alpha.*resultGUI_matRad.physicalDose)), hold on, plot(matRad_calcIDD(sqrt(resultGUI_matRad.beta).*resultGUI_matRad.physicalDose))
% plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.alpha.*resultGUI_TOPAS_LEM1.physicalDose))
% plot(matRad_calcIDD(sqrt(resultGUI_TOPAS_LEM1.beta).*resultGUI_TOPAS_LEM1.physicalDose))
% legend({'matRad alpha*dose','matRad sqrt(beta)*dose','TOPAS LEM1 alpha*dose','TOPAS LEM1 sqrt(beta)*dose'})
% %%
% figure, plot(matRad_calcIDD(resultGUI_matRad.alpha)), hold on, plot(matRad_calcIDD(resultGUI_matRad.beta))
% plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.alpha))
% plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.beta))
% legend({'matRad alpha','matRad beta','TOPAS LEM1 alpha','TOPAS LEM1 beta'})
% %%
% 
% figure, plot(matRad_calcIDD(resultGUI_matRad.alpha,'y',1)), hold on, plot(matRad_calcIDD(resultGUI_matRad.beta,'y',1))
% plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.alpha,'y',1))
% plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.beta,'y',1))
% legend({'matRad alpha','matRad beta','TOPAS LEM1 alpha','TOPAS LEM1 beta'})
% 
% figure, plot(matRad_calcIDD(resultGUI_matRad.alpha.*resultGUI_matRad.physicalDose,'y',1)), hold on, plot(matRad_calcIDD(sqrt(resultGUI_matRad.beta).*resultGUI_matRad.physicalDose,'y',1))
% plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.alpha.*resultGUI_TOPAS_LEM1.physicalDose,'y',1))
% plot(matRad_calcIDD(sqrt(resultGUI_TOPAS_LEM1.beta).*resultGUI_TOPAS_LEM1.physicalDose,'y',1))
% legend({'matRad alpha*dose','matRad sqrt(beta)*dose','TOPAS LEM1 alpha*dose','TOPAS LEM1 sqrt(beta)*dose'})
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% figure, hold on,
% plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.physicalDose),'--')
% plot(matRad_calcIDD(resultGUI_TOPAS_LEM1.RBExD))
% plot(matRad_calcIDD(resultGUI_TOPAS_direct_old.physicalDose),'--')
% plot(matRad_calcIDD(resultGUI_TOPAS_direct_old.RBExD))
% legend({'TOPAS LEM1 physDose','TOPAS LEM1','TOPAS amtrack physDose','TOPAS amtrack'})
% %%
% 
% figure, plot(matRad_calcIDD(resultGUI_TOPAS_direct.physicalDose)), hold on, plot(matRad_calcIDD(resultGUI_TOPAS_fit.physicalDose))
% xlim([40 100])
% legend({'matRad','TOPAS'})
% % RBE = readBinData('A:\matRad_VARBEmerge\topas\MCrun\score_matRad_plan_field1_run1_RBE.bin',ct.cubeDim);
% % physDose = readBinData('A:\matRad_VARBEmerge\topas\MCrun\score_matRad_plan_field1_run1_physicalDose.bin',ct.cubeDim);
% % 
% % figure, plot(matRad_calcIDD(physDose,'y')), hold on, plot(matRad_calcIDD(physDose.*RBE,'y'))
% 
% %%
% % resultGUI_TOPAS_RBE = matRad_calcDoseDirectMC(ct,stf,pln,cst,weights,1e5);
% %%
% % dijTOPAS = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst);
% % resultGUI_TOPAS = matRad_fluenceOptimization(dijTOPAS,cst,pln);
% 
% figure, plot(matRad_calcIDD(resultGUI_matRad.physicalDose)), hold on,...
% plot(matRad_calcIDD(resultGUI_TOPAS.physicalDose))
% if strcmp(pln.radiationMode,'protons')
%     plot(matRad_calcIDD(resultGUI_MCsquare.physicalDose)),...
% end
% legend({'matRad','TOPAS','MCsquare'},'Location','northwest')
% xlim([10 100])
% ylabel('carbon physicalDose')
% %%
% figure,plot(matRad_calcIDD(resultGUI_matRad.RBExD)) , hold on,...
% plot(matRad_calcIDD(resultGUI_TOPAS.RBE.*resultGUI_TOPAS.physicalDose))
% legend({'matRad','TOPAS'},'Location','northwest')
% xlim([10 100])
% ylabel('carbon RBExD')
% %%
% loyolagreen = 1/255*[0,104,87];
% blue = [0, 0.4470, 0.7410];
% orange = [0.8500, 0.3250, 0.0980];
% 
% figure, plot(matRad_calcIDD(resultGUI_matRad.physicalDose),'--','Color',blue), hold on,...
% plot(matRad_calcIDD(resultGUI_matRad.RBExD),'Color',blue)
% plot(matRad_calcIDD(resultGUI_TOPAS.physicalDose),'--','Color',orange)
% plot(matRad_calcIDD(resultGUI_TOPAS.RBExD),'Color',orange)
% xlim([10 90])
% ylim([0 0.11])
% legend({'matRad physDose','matRad RBExD','TOPAS physDose','TOPAS RBExD'},'Location','northwest')
% %% FITS
% %% carbon
% data = importdata('CarbonWater_mgcm-2.txt');
% range = [39 48];
% Range = data.data(range(1):range(2),5)/100;
% Energy = data.data(range(1):range(2),1);
% 
% meanEnergy = @(x) 10.55 * x^0.6409 + 14.14;
% 
% %% helium
% data = importdata('HeliumWater_mgcm-2.txt');
% range = [40 45];
% Range = data.data(range(1):range(2),5)*100;
% Energy = data.data(range(1):range(2),1);
% 
% meanEnergy = @(x) 9.466* x.^0.5615 + 1.719;
% meanEnergy_old = @(x) 10.12* x.^0.5519 - 0.5487;




