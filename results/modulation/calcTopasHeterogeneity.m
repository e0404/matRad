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
stf = matRad_generateStfPencilBeam(pln,ct);
% stf = matRad_generateStf(ct,cst,pln);

weights = ones(1,sum([stf(:).totalNumOfBixels]));
resultGUI_matRad = matRad_calcDoseDirect(ct,stf,pln,cst,weights);

%% Analytical heterogeneity correction
cstHetero = matRad_cstHeteroAutoassign(cst);
pln.propHeterogeneity.calcHetero = true;
resultGUI_matRad_hetero = matRad_calcDoseDirect(ct,stf,pln,cstHetero,resultGUI_matRad.w);

%% TOPAS heterogeneity correction
pln.propMC.proton_engine = 'TOPAS';
resultGUI_TOPAS = matRad_calcDoseDirectMC(ct,stf,pln,cst,weights,1e5);
%%
% dijTOPAS = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst);
% resultGUI_TOPAS = matRad_fluenceOptimization(dijTOPAS,cst,pln);


%%
samples = [5, 10, 50];

%%% matRad
numOfSamples = max(samples);
for i = 1:numOfSamples
    ct_mod = matRad_modulateDensity(ct,cst,800);
    i
    resultGUI_mod{i} = matRad_calcDoseDirect(ct_mod,stf,pln,cst,resultGUI_matRad.w);                                                                      
end
%%
%%% TOPAS
pln.propMC.proton_engine = 'TOPAS';
for i = 1:numOfSamples
    ct_mod = matRad_modulateDensity(ct,cst,800);
    i
    resultGUI_modTOPAS{i} = matRad_calcDoseDirectMC(ct_mod,stf,pln,cst,weights,1e5);                                                                     
end
%%
%%% TOPAS aufgeteilt
for s = samples
    for i = 1:s
        ct_mod = matRad_modulateDensity(ct,cst,800);
        i
        resultGUI_modTOPASaufgeteilt.(['samples',num2str(s)]){i} = matRad_calcDoseDirectMC(ct_mod,stf,pln,cst,weights,1e5/s);
    end
end
%%

%%% matRad
for s = samples
    resultGUI_matRad_mod.(['physicalDose',num2str(s)]) = zeros(160,160,160);
    for i = 1:s
        resultGUI_matRad_mod.(['physicalDose',num2str(s)]) = resultGUI_matRad_mod.(['physicalDose',num2str(s)]) + resultGUI_mod{i}.physicalDose/s;
    end
end
%%
%%% TOPAS
for s = samples
    resultGUI_TOPAS_mod.(['physicalDose',num2str(s)]) = zeros(160,160,160);
    for i = 1:s
        resultGUI_TOPAS_mod.(['physicalDose',num2str(s)]) = resultGUI_TOPAS_mod.(['physicalDose',num2str(s)]) + resultGUI_modTOPAS{i}.physicalDose/s;
    end
end

%%% TOPAS aufgeteilt
for s = samples
    resultGUI_TOPAS_mod_aufgeteilt.(['physicalDose',num2str(s)]) = zeros(160,160,160);
    for i = 1:s
        resultGUI_TOPAS_mod_aufgeteilt.(['physicalDose',num2str(s)]) = resultGUI_TOPAS_mod_aufgeteilt.(['physicalDose',num2str(s)]) + resultGUI_modTOPASaufgeteilt.(['samples',num2str(s)]){i}.physicalDose/s;
    end
end

%%
% save homogeneousPhantom.mat -v7.3
%%
figure, hold on, ...
    plot(matRad_calcIDD(resultGUI_matRad.physicalDose)), ...
    plot(matRad_calcIDD(resultGUI_TOPAS.physicalDose)),...
    
    plot(matRad_calcIDD(resultGUI_matRad_hetero.physicalDose)),...
    
    plot(matRad_calcIDD(resultGUI_matRad_mod.physicalDose50)),...

    plot(matRad_calcIDD(resultGUI_TOPAS_mod.physicalDose50),':'),...
    plot(matRad_calcIDD(resultGUI_TOPAS_mod_aufgeteilt.physicalDose50),'--')

legend({'matRad','TOPAS','matRad analytical','matRad modulated','TOPAS modulated','TOPAS modulated split histories'},'Location','northwest')
xlim([10 100])


figure, hold on, ...
    plot(matRad_calcIDD(resultGUI_matRad_mod.physicalDose5)),...
    plot(matRad_calcIDD(resultGUI_matRad_mod.physicalDose10)),...
    plot(matRad_calcIDD(resultGUI_matRad_mod.physicalDose50)),...

    plot(matRad_calcIDD(resultGUI_TOPAS_mod.physicalDose5),':'),...
    plot(matRad_calcIDD(resultGUI_TOPAS_mod.physicalDose10),':'),...
    plot(matRad_calcIDD(resultGUI_TOPAS_mod.physicalDose50),':'),...

    plot(matRad_calcIDD(resultGUI_TOPAS_mod_aufgeteilt.physicalDose5),'--'),...
    plot(matRad_calcIDD(resultGUI_TOPAS_mod_aufgeteilt.physicalDose10),'--'),...
    plot(matRad_calcIDD(resultGUI_TOPAS_mod_aufgeteilt.physicalDose50),'--')

legend({'matRad 5 samples','matRad 10 samples','matRad 50 samples','TOPAS 5 samples','TOPAS 10 samples','TOPAS 50 samples','TOPAS 5 samples split histories','TOPAS 10 samples split histories','TOPAS 50 samples split histories'},'Location','northwest')
xlim([10 100])

% figure, plot(matRad_calcIDD(resultGUI_matRad.physicalDose)), hold on, plot(matRad_calcIDD(resultGUI_matRad_hetero.physicalDose))
% for s = samples
%     plot(matRad_calcIDD(resultGUI_matRad_mod.(['physicalDose',num2str(s)])))
% end
% 
% title('Homogeneous Phantom')
% legend({'matRad','matRad analytisch korrigiert','matRad modulated 100 samples','matRad modulated 10 samples'},'Location','northwest')
% xlim([50 90])
%%
% load CalcHomogeneousLungPhantom

% %%
% load CalcHeterogeneousLungPhantom
% figure, plot(matRad_calcIDD(resultGUI_matRad.physicalDose)), hold on, plot(matRad_calcIDD(resultGUI_matRad_hetero.physicalDose))
% plot(matRad_calcIDD(resultGUI_TOPAS.physicalDose))
% plot(matRad_calcIDD(resultGUI_MC_mod.physicalDose))
% title('Heterogeneous Phantom')
% legend({'matRad','matRad analytisch korrigiert','TOPAS','matRad modulated'},'Location','northwest')
% xlim([50 90])


