clear

%% Patient Data Import
load('BOXPHANTOM_LUNG_LARGE');
% cstDev = cst;

% for i = 1:3
%     cst{i,6} = struct;
%     cst{i,6}.dose  = cstDev{i,6}{1,1}.parameters{1,1};
%     cst{i,6}.penalty = cstDev{i,6}{1,1}.penalty;
%     cst{i,6}.robustness = 'none';
% end
% 
% cst{1,6}.type = 'square overdosing';
% cst{2,6}.type = 'square deviation';
% cst{3,6}.type = 'square overdosing';
% clear cstDev

%% Treatment Plan

pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.machine         = 'generic_TOPAS_cropped_APM';

%%
modelName           = 'none';
quantityOpt         = 'PhysicalDose';

pln.propOpt.bioOptimization = 'none';

%%
pln.numOfFractions = 1;

pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 10;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen'); % optimize on the nominal scenario

pln.heterogeneity.calcHetero = true;
pln.heterogeneity.type = 'complete'; % optional
pln.heterogeneity.useDoseCurves = true;
pln.heterogeneity.useOrgDepths = false;

stf = matRad_generateStf(ct,cst,pln,param);
%%

stf1 = stf;

%% Selecting the center ray
for r = 1:length(stf1.ray)
    if stf1.ray(r).rayPos_bev(1) == 0 && stf1.ray(r).rayPos_bev(2) == 0 && stf1.ray(r).rayPos_bev(3) == 0
        break
    end
end

stf.ray = stf1.ray(r);
energyIx = round(numel(stf1.ray(r).energy)/2);
stf.ray.energy = stf1.ray(r).energy(energyIx);
stf.ray.focusIx = stf1.ray(r).focusIx(energyIx);
stf.ray.rangeShifter = stf1.ray(r).rangeShifter(energyIx);

stf.numOfRays = 1;
stf.numOfBixelsPerRay = 1;
stf.totalNumOfBixels = 1;

%%
% addpath(genpath('A:\matRad_VARRBE'))
% rmpath(genpath('A:\matRad_Noa'))
% 
% pln.propMC.proton_engine = 'TOPAS';
% pln.machine         = 'generic_TOPAS_cropped';
% numOfSamples = 50;
% 
% dij = matRad_calcParticleDose(ct,stf,pln,cst,param);
% resultGUI = matRad_fluenceOptimization(dij,cst,pln,param);
% 
% rmpath(genpath('A:\matRad_VARRBE'))
% addpath(genpath('A:\matRad_Noa'))
% 
% for i = 1:numOfSamples
%     ct_mod = matRad_modulateDensity(ct,cst,800);
%     resultGUI_mod{i} = matRad_calcDoseDirect(ct_mod,stf,pln,cst,resultGUI.w);
% end
% 
% resultGUI_matRad_mod.physicalDose = zeros(160,160,160);
% for i = 1:numOfSamples
%     resultGUI_matRad_mod.physicalDose = resultGUI_matRad_mod.physicalDose + resultGUI_mod{i}.physicalDose/numOfSamples;
% end
% 
% 
% %%
% matRad_compareDose(resultGUI.physicalDose, resultGUI_matRad_mod.physicalDose, ct, cst, [0, 1, 0] , 'off', pln, [2, 2], 1, 'global');





%% Analytical heterogeneity correction

% dij = matRad_calcParticleDose(ct,stf,pln,cst,param);
% resultGUI = matRad_fluenceOptimization(dij,cst,pln,param);
resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,1,param);

cstHetero = matRad_cstHeteroAutoassign(cst);

resultGUI_hetero = matRad_calcDoseDirect(ct,stf,pln,cstHetero,resultGUI.w,param);



%% TOPAS heterogeneity correction

rmpath(genpath('A:\matRad_VARRBE'))
addpath(genpath('A:\matRad_Noa'))

pln.propMC.proton_engine = 'TOPAS';
pln.machine         = 'generic_TOPAS_cropped';
numOfSamples = 50;

resultGUI_MC = matRad_calcDoseDirectMC(ct,stf,pln,cst,resultGUI.w,1e5);

for i = 1:numOfSamples
    ct_mod = matRad_modulateDensity(ct,cst,800);
    resultGUI_mod{i} = matRad_calcDoseDirectMC(ct_mod,stf,pln,cst,resultGUI.w,1e5);                                                                      
end

resultGUI_MC_mod.physicalDose = zeros(160,160,160);
for i = 1:numOfSamples
    resultGUI_MC_mod.physicalDose = resultGUI_MC_mod.physicalDose + resultGUI_mod{i}.physicalDose/numOfSamples;
end

%save('TopasHeterogeneity')
%% Compare Dose
% pause(10);
%system('shutdown -s -t 3600')
% 
% matRad_compareDose(resultGUI_MC.physicalDose, resultGUI_MC_mod.physicalDose, ct, cst, [0, 1, 0] , 'off', pln, [2, 2], 1, 'global');
% 
% %%
% 
% 
% 
% 
% %% Create large Lung
% lungBox = 20:140;
% for i = lungBox
%     for j = lungBox
%         for k = lungBox
%             ct.cube{1}(i,j,k) = 0.4;
%         end
%     end
% end
% 
% body = find(ct.cube{1} == 0.4);
% box = 50:110;
% for i = box
%     for j = box
%         for k = box
%             ct.cube{1}(i,j,k) = 1;
%         end
%     end
% end
% lung = find(ct.cube{1} == 0.4);
% ct.cubeHU{1}(:) = 1024*(ct.cube{1}(:)-1);
% 
% cst{1,4}{1,1} = body;
% cst{3,4}{1,1} = lung;
%%

load('S003_3.mat')

ct_mod = matRad_modulateDensity(ct,cst,800);
lungSelect = strcmp(cst(:,2),'Lunge bds');
tumorSelect = strcmp(cst(:,2),'PTV');

selection = lungSelect + tumorSelect;

cst{lungSelect,5}.visibleColor = [0 0 1];
cst{tumorSelect,5}.visibleColor = [1 0 0];

% ct.cubeHU2{1} = ct.cubeHU{1}
f = figure;
set(f, 'Units', 'normalized', 'Position', [0.3, 0.3, 0.5, 0.3]);

for i = 104

slice = i;
subplot(1,2,1)
matRad_plotCtSlice(gca,ct.cubeHU,1,1,slice)
matRad_plotVoiContourSlice(gca,cst,ct.cubeHU,1,selection,1,slice);
xlabel(sprintf('y [voxel / %i mm]',ct.resolution.x))
ylabel(sprintf('x [voxel / %i mm]',ct.resolution.y))
xlim([40 120])
ylim([20 140])

title('Original HU cube')
plotCustomLegend
set(gca, 'XAxisLocation', 'top')
camroll(90)

subplot(1,2,2)
matRad_plotCtSlice(gca,ct_mod.cubeHU,1,1,slice)
matRad_plotVoiContourSlice(gca,cst,ct_mod.cubeHU,1,selection,1,slice);
xlabel(sprintf('y [voxel / %i mm]',ct.resolution.x))
ylabel(sprintf('x [voxel / %i mm]',ct.resolution.y))
xlim([40 120])
ylim([20 140])

title('HU cube with sampled lung')
plotCustomLegend
set(gca, 'XAxisLocation', 'top')
camroll(90)

%waitforbuttonpress;

end





%%
clear
load('TopasHeterogeneity.mat')
%%
figure
homogeneous = matRad_calcIDD(resultGUI.physicalDose,'y');
normalization = 1/max(homogeneous);

heterogeneousMatRad = matRad_calcIDD(resultGUI_hetero.physicalDose,'y');
homogeneousMC = matRad_calcIDD(resultGUI_MC.physicalDose,'y');
heterogeneousMC = matRad_calcIDD(resultGUI_MC_mod.physicalDose,'y');

linewidth = 1.5;
plot(homogeneous*normalization,'Linewidth',linewidth);
hold on
plot(heterogeneousMatRad*normalization,'Linewidth',linewidth);
plot(homogeneousMC*normalization,'Linewidth',linewidth);
plot(heterogeneousMC*normalization,'Linewidth',linewidth);

patch([20 20 50 50],[-1 1.5 1.5 -1],'black','facealpha',0.2)

xlabel(sprintf('y [voxel / %i mm]',ct.resolution.y))
ylabel('Dose [normalized to matRad]')
legend({'matRad homogeneous','matRad corrected','TOPAS homogeneous','TOPAS corrected','Lung'},'Location','northwest')

y1 = 0.35;
y2 = 1.15;
x1 = 70;
x2 = 82;
patch([x1 x1 x2 x2],[y1 y2 y2 y1],'black','facealpha',0,'HandleVisibility','off')

xlim([20 100])
ylim([0 1.2])

ax1 = gca;
f2 = figure();
set(f2,'Position',[1575 918 275 420])
ax2 = copyobj(ax1,f2);
xlim([x1 x2])
ylim([y1 y2])



