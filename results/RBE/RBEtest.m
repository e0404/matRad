% clear

% Patient Data Import
load('BOXPHANTOM.mat');
ct.cube{1} = ones(160,160,160);
ct.cubeHU{1} = zeros(160,160,160);
cst{1,4}{1,1} = [1:prod([160,160,160])]';
% Treatment Plan

pln.radiationMode   = 'helium';     % either photons / protons / carbon
% pln.machine         = 'HITgenericRIFI3MMTOPAS_Water_forNiklasAPM';
pln.machine         = 'Generic';

modelName           = 'none';
% modelName           = 'MCN';
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

stf.ray.energy = machine.data(60).energy;
resultGUI_matRad_3 = matRad_calcDoseDirect(ct,stf,pln,cst,weights);

%% TOPAS heterogeneity correction
if strcmp(pln.radiationMode,'protons')
    pln.propMC.proton_engine = 'MCsquare';
    resultGUI_MCsquare = matRad_calcDoseDirectMC(ct,stf,pln,cst,1,1e5);
end
%%
% pln.propMC.carbon_engine = 'TOPAS';
tic
resultGUI_TOPAS_direct = matRad_calcDoseDirectMC(ct,stf,pln,cst,1,1e3);
toc
%%
figure, plot(matRad_calcIDD(resultGUI_matRad_3.physicalDose)), hold on, plot(matRad_calcIDD(resultGUI_TOPAS_direct.physicalDose))
%%


figure, plot(matRad_calcIDD(resultGUI_TOPAS_direct.physicalDose)), hold on, plot(matRad_calcIDD(resultGUI_TOPAS_fit.physicalDose))
xlim([40 100])
legend({'matRad','TOPAS'})
% RBE = readBinData('A:\matRad_VARBEmerge\topas\MCrun\score_matRad_plan_field1_run1_RBE.bin',ct.cubeDim);
% physDose = readBinData('A:\matRad_VARBEmerge\topas\MCrun\score_matRad_plan_field1_run1_physicalDose.bin',ct.cubeDim);
% 
% figure, plot(matRad_calcIDD(physDose,'y')), hold on, plot(matRad_calcIDD(physDose.*RBE,'y'))

%%
% resultGUI_TOPAS_RBE = matRad_calcDoseDirectMC(ct,stf,pln,cst,weights,1e5);
%%
% dijTOPAS = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst);
% resultGUI_TOPAS = matRad_fluenceOptimization(dijTOPAS,cst,pln);

figure, plot(matRad_calcIDD(resultGUI_matRad.physicalDose)), hold on,...
plot(matRad_calcIDD(resultGUI_TOPAS.physicalDose))
if strcmp(pln.radiationMode,'protons')
    plot(matRad_calcIDD(resultGUI_MCsquare.physicalDose)),...
end
legend({'matRad','TOPAS','MCsquare'},'Location','northwest')
xlim([10 100])
ylabel('carbon physicalDose')
%%
figure,plot(matRad_calcIDD(resultGUI_matRad.RBExD)) , hold on,...
plot(matRad_calcIDD(resultGUI_TOPAS.RBE.*resultGUI_TOPAS.physicalDose))
legend({'matRad','TOPAS'},'Location','northwest')
xlim([10 100])
ylabel('carbon RBExD')
%%
loyolagreen = 1/255*[0,104,87];
blue = [0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];

figure, plot(matRad_calcIDD(resultGUI_matRad.physicalDose),'--','Color',blue), hold on,...
plot(matRad_calcIDD(resultGUI_matRad.RBExD),'Color',blue)
plot(matRad_calcIDD(resultGUI_TOPAS.physicalDose),'--','Color',orange)
plot(matRad_calcIDD(resultGUI_TOPAS.RBExD),'Color',orange)
xlim([10 90])
ylim([0 0.11])
legend({'matRad physDose','matRad RBExD','TOPAS physDose','TOPAS RBExD'},'Location','northwest')
%% FITS
%% carbon
data = importdata('CarbonWater_mgcm-2.txt');
range = [39 48];
Range = data.data(range(1):range(2),5)/100;
Energy = data.data(range(1):range(2),1);

meanEnergy = @(x) 10.55 * x^0.6409 + 14.14;

%% helium
data = importdata('HeliumWater_mgcm-2.txt');
range = [40 45];
Range = data.data(range(1):range(2),5)*100;
Energy = data.data(range(1):range(2),1);

meanEnergy = @(x) 9.466* x.^0.5615 + 1.719;
meanEnergy_old = @(x) 10.12* x.^0.5519 - 0.5487;




