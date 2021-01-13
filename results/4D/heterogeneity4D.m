%% 4D dose calculation workflow
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In this example we will show 
% (i) the structure of 4D data within matRad
% (ii) how to perform standard treatment planning
% (iii) how to run a dose recalculation considering interplay effects 
 
%% set matRad runtime configuration
clear
matRad_rc

%% Load data, add generic 4D information, and display 'moving' geometry
% load BOXPHANTOM_LUNG_LARGE_2e-1.mat
load BOXPHANTOM.mat

%%

amplitude    = [0 5 0]; % [voxels]
numOfCtScen  = 5;
motionPeriod = 2.5; % [s] 

[ct,cst] = matRad_addMovement(ct, cst,motionPeriod, numOfCtScen, amplitude);
% Set up a plan, compute dose influence on all phases, conventional optimization
% meta information for treatment plan

%%
% figure
% for i = 1:length(ct.cube)
%    hold on
%    imagesc(ct.cube{i}(:,:,80));
%    waitforbuttonpress;
% end
%%


pln.numOfFractions  = 15;
pln.radiationMode   = 'carbon';           % either photons / protons / helium / carbon
pln.machine = 'HITgenericRIFI3MMTOPAS_Water_forNiklas';
% pln.radiationMode   = 'protons';  
% pln.machine         = 'generic_TOPAS_cropped_APM';


% beam geometry settings
pln.propStf.bixelWidth      = 7; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longitudinalSpotSpacing = 10;      % only relevant for HIT machine, not generic
pln.propStf.gantryAngles    = 90; 
pln.propStf.couchAngles     = 0; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

%optimization settings
pln.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

% quantityOpt  = 'RBExD';     % options: physicalDose, effect, RBExD
% modelName    = 'MCN';             % none: for photons, protons, carbon            % constRBE: constant RBE 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions
quantityOpt  = 'RBExD';
modelName    = 'LEM';
scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen' 

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,scenGenType);

%%
% generate steering file
stf = matRad_generateStf(ct,cst,pln);
% stf = matRad_generateStfPencilBeam(pln,ct);

%% 
% dose calculation
% weights = ones(1,sum([stf(:).totalNumOfBixels]));
% resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,weights);
dij = matRad_calcParticleDose(ct,stf,pln,cst);

%% 
% inverse planning for imrt on a static CT
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% 
% post processing
% This step is necessary to remove beam spots with too few particles that
% cannot not be delivered, dose is recalculated accordingly
resultGUI = matRad_postprocessing(resultGUI, dij, pln, cst, stf) ; 

%%
% pln.propHeterogeneity.calcHetero = true;
% cstHetero = matRad_cstHeteroAutoassign(cst);
%%
% dijHetero = matRad_calcParticleDose(ct,stf,pln,cstHetero);
% resultGUI_hetero = matRad_fluenceOptimization(dijHetero,cstHetero,pln);
% resultGUI_hetero = matRad_postprocessing(resultGUI_hetero, dijHetero, pln, cstHetero, stf) ; 
%%
% resultGUI_heteroDirect = matRad_calcDoseDirect(ct,stf,pln,cstHetero,resultGUI.w);
%%
% matRad_compareDose(resultGUI.physicalDose, resultGUI_hetero.physicalDose , ct, cst, [1, 1, 0] , 'off', pln, [2, 2], 1, 'global');
% figure, plot(matRad_calcIDD(resultGUI.physicalDose,'x')), hold on, plot(matRad_calcIDD(resultGUI_hetero.physicalDose,'x')), plot(matRad_calcIDD(resultGUI_heteroDirect.physicalDose,'x'),'--'), legend({'homogeneous','with correction dij','with correction direct'})
% xlim([60 150])
% figure, plot(matRad_calcIDD(resultGUI_hetero.physicalDose,'x')), hold on, plot(matRad_calcIDD(resultGUI_heteroDirect.physicalDose,'x')), legend({'with correction dij','with correction direct'})
% figure, plot(matRad_calcIDD(resultGUI_heteroDirect.physicalDose,'x'))
%%
% calc 4D dose
% make sure that the correct pln, dij and stf are loeaded in the workspace
[resultGUI, timeSequence] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI); 

%%
% TOPAS calculation
pln.propMC.proton_engine = 'TOPAS';
resultGUI_topas = matRad_calcDoseDirectMC(ct,stf,pln,cst,timeSequence.phaseMatrix,1e3);
 
%%
save 4DprotonTOPAScroppedRBE
figure, plot(matRad_calcIDD(resultGUI.phaseDose{1})), hold on, plot(matRad_calcIDD(resultGUI_topas.phaseDose{1}))

%%
figure, hold on,
for i = 1:5
    plot(matRad_calcIDD(resultGUI.phaseDose{i},'y'))
%     imagesc(resultGUI.phaseDose{i}(:,:,80))

%     imagesc(ct.cube{1,i}(:,:,80))
    waitforbuttonpress;
end

%%
figure, hold on,
for i = 1:5
    plot(matRad_calcIDD(resultGUI_topas.phaseDose{i},'y'))
%     imagesc(resultGUI.phaseDose{i}(:,:,80))

%     imagesc(ct.cube{1,i}(:,:,80))
    waitforbuttonpress;
end

%%
figure, plot(matRad_calcIDD(resultGUI.physicalDose,'y')), hold on, plot(matRad_calcIDD(resultGUI_topas.physicalDose,'y'))



%%
%plot the result in comparison to the static dose
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z); 

figure 

subplot(2,2,1)
imagesc(resultGUI.physicalDose(:,:,slice)),colorbar, colormap(jet);
title('static dose distribution [Gy (RBE)]')
axis equal

subplot(2,2,3)
imagesc(resultGUI.accPhysicalDose(:,:,slice)),colorbar, colormap(jet); 
title('accumulated (4D) dose distribution [Gy (RBE)]')
axis equal

subplot(2,2,2)
imagesc(resultGUI.physicalDose(:,:,slice) - resultGUI.accPhysicalDose(:,:,slice)) ,colorbar, colormap(jet); 
title('static dose distribution - accumulated (4D) dose distribution [Gy (RBE)]')

axis equal
% 
