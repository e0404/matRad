% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% might need a circular phantom
% load BOXPHANTOM.mat
%load LIVER_HP.mat
% load LIVERSingleC12.mat
load LiverProt.mat

% meta information for treatment plan


% photon bit pln(1)
pln(1).bixelWidth      = 3; % [mm] / also corresponds to lateral spot spacing for particles
pln(1).gantryAngles    = [0:72:359]; % [?]
pln(1).couchAngles     = zeros( numel(pln(1).gantryAngles),1); % [?]
pln(1).numOfBeams      = numel(pln(1).gantryAngles);
pln(1).numOfVoxels     = prod(ct.cubeDim);
pln(1).isoCenter       = ones(pln(1).numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln(1).voxelDimensions = ct.cubeDim;
pln(1).radiationMode   = 'photons';     % either photons / protons / carbon

pln(1).bioOptimization = 'test_effect';         % none_physicalDose: physical optimization;                              constRBE_RBExD; constant RBE of 1;  
                                             % MCN_effect; McNamara-variable RBE model for protons (effect based)     MCN_RBExD; McNamara-variable RBE model for protons (RBExD) based
                                             % WED_effect; Wedenberg-variable RBE model for protons (effect based)    WED_RBExD; Wedenberg-variable RBE model for protons (RBExD) based
                                             % LEM_effect: effect-based optimization;                                 LEM_RBExD: optimization of RBE-weighted dose
                                             % test_effect: effect based optimization for photons                     test  
pln(1).scenGenType     = 'nomScen';             % scenario creation type'nomScen'  'wcScen' 'impScen' 'rndScen'                                              
pln(1).numOfFractions  = 30;
pln(1).runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(1).runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(1).machine         = 'Generic';
pln(1).robOpt          = false;
quantityOpt         = 'physicalDose';     % options: physicalDose, constRBE, effect, RBExD
modelName           = 'none';             % none: for photons, protons, carbon                                    MCN: McNamara-variable RBE model for protons
                                          % WED: Wedenberg-variable RBE model for protons                         LEM: Local Effect Model for carbon ions
% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% set plan uncertainties for robust optimization
pln(1).multScen = matRad_multScen(ct,pln(1).scenGenType); 


% proton pln(2)
pln(2).bixelWidth      = 3; % [mm] / also corresponds to lateral spot spacing for particles
pln(2).gantryAngles    = [300 ]; % [?]
pln(2).couchAngles     = [ 0]; % [?]
pln(2).numOfBeams      = numel(pln(2).gantryAngles);
pln(2).numOfVoxels     = prod(ct.cubeDim);
pln(2).isoCenter       = ones(pln(2).numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln(2).voxelDimensions = ct.cubeDim;
pln(2).radiationMode   = 'protons';     % either photons / protons / carbon

pln(2).bioOptimization = 'MCN_effect';        % none_physicalDose: physical optimization;                              constRBE_RBExD; constant RBE of 1.1;  
                                             % MCN_effect; McNamara-variable RBE model for protons (effect based)     MCN_RBExD; McNamara-variable RBE model for protons (RBExD) based
                                             % WED_effect; Wedenberg-variable RBE model for protons (effect based)    WED_RBExD; Wedenberg-variable RBE model for protons (RBExD) based
                                             % LEM_effect: effect-based optimization;                                 LEM_RBExD: optimization of RBE-weighted dose
pln(2).scenGenType     = 'nomScen';             % scenario creation type'nomScen'  'wcScen' 'impScen' 'rndScen'                                                          
pln(2).numOfFractions  = 30;
pln(2).runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(2).runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(2).machine         = 'Generic';
pln(2).robOpt          = false;
quantityOpt         = 'physicalDose';     % options: physicalDose, constRBE, effect, RBExD
modelName           = 'none';             % none: for photons, protons, carbon                                    MCN: McNamara-variable RBE model for protons
                                          % WED: Wedenberg-variable RBE model for protons                         LEM: Local Effect Model for carbon ions
% retrieve bio model parameters
pln(2).bioParam = matRad_bioModel(pln(2).radiationMode,quantityOpt, modelName);

% set plan uncertainties for robust optimization
pln(2).multScen = matRad_multScen(ct,pln(2).scenGenType); 

% carbon bit pln(1)
pln(3).bixelWidth      = 3; % [mm] / also corresponds to lateral spot spacing for particles
pln(3).gantryAngles    = [300 ]; % [?]
pln(3).couchAngles     = [0 ]; % [?]
pln(3).numOfBeams      = numel(pln(3).gantryAngles);
pln(3).numOfVoxels     = prod(ct.cubeDim);
pln(3).isoCenter       = ones(pln(3).numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln(3).voxelDimensions = ct.cubeDim;
pln(3).radiationMode   = 'carbon';     % either photons / protons / carbon
%pln.radiationMode   = {'carbon','carbon','photons'};     % either photons / protons / carbon
pln(3).bioOptimization = 'LEM_effect';        % none_physicalDose: physical optimization;                              constRBE_RBExD; constant RBE of 1.1;  
                                             % MCN_effect; McNamara-variable RBE model for protons (effect based)     MCN_RBExD; McNamara-variable RBE model for protons (RBExD) based
                                             % WED_effect; Wedenberg-variable RBE model for protons (effect based)    WED_RBExD; Wedenberg-variable RBE model for protons (RBExD) based
                                             % LEM_effect: effect-based optimization;                                 LEM_RBExD: optimization of RBE-weighted dose
pln(3).scenGenType     = 'nomScen';             % scenario creation type'nomScen'  'wcScen' 'impScen' 'rndScen'                                                          
pln(3).numOfFractions  = 1;
pln(3).runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(3).runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(3).machine         = 'Generic';
pln(3).robOpt          = false;
quantityOpt         = 'physicalDose';     % options: physicalDose, constRBE, effect, RBExD
modelName           = 'none';             % none: for photons, protons, carbon                                    MCN: McNamara-variable RBE model for protons
                                          % WED: Wedenberg-variable RBE model for protons                         LEM: Local Effect Model for carbon ions
% retrieve bio model parameters
pln(3).bioParam = matRad_bioModel(pln(3).radiationMode,quantityOpt, modelName);

% set plan uncertainties for robust optimization
% retrieve scenarios for dose calculation and optimziation
pln(3).multScen = matRad_multScen(ct,pln(3).scenGenType); 




%% generate stfs

% stf = matRad_generateStf(ct,cst,pln);
% param.mixedMod = [1 2 4] if you wanna choose 1st,2nd 4th plns 

param.mixedMod = [1 2];

stf = matRad_generateStf(ct,cst,pln,param); 

%% calc multiple modality dose
% here generate Superpln for the Optimizer, newPln for modded meta for
% modalities, as its adjusted to "max complexity"
pln1= pln;
[dij, dijSepa, superPln] = matRad_calcCombiDose(ct,stf,pln,cst,param);
%% preconditioning
% dij_copy1 = dij;
% [dij,arbit] = matRad_mixedModPreconditioner(dij, stf, param); 
% dij.wConditionCorrection = arbit;
% review the need for it 
%%  FRACTIONATION SETUP IN THIS BRANCH 
% linear scaling of weights in back projection 
% also change the objectives from fraction based to total ("total effect")
% after total effect split into smaller sets 
%

% 
%% fluence opti uses some parts of pln for options change basic pipeline for it
resultGUI = matRad_fluenceOptimization(dij,cst,superPln,param, dijSepa);
result_copy  = resultGUI;
%% keep a copy 

dij_copy  = dij;
pln_copy = pln;
result_copy = resultGUI;
%%
% 
% %cube =  result_copy{1}.physicalDose + result_copy{2}.physicalDose ;
%  resultGUI = result_copy{1};
%  dij = dijSepa{1};
%  pln = pln(1);
%  %resultGUI.physicalDose_beam1 = result_copy{3}.RBExD;
%  resultGUI.physicalDose_beam2 = result_copy{2}.RBExD;
%  resultGUI.physicalDose = cube;
%%
cube =  result_copy{1}.RBExD + result_copy{2}.RBExD;% + result_copy{3}.RBExD ;
i=2;
resultGUI =result_copy{i};
dij = dijSepa{i};
pln = pln_copy(i);

 resultGUI.physicalDose = cube;
%%
matRadGUI
% 
% % 
% 
%  %%
% % wInit = ones(sum([stf.totalNumOfBixels]),1);
% % wInit = resultGUI.w;
% % 
% % wInitxRay = wInit;
% % wInitProt = wInit;
% % wInitCarb = wInit;
% % wTot      = wInit;
% % off=1;
% % for i=1:numel(stf)
% %     if strcmp(stf(i).radiationMode, 'photons')
% %         off=off+stf(i).numOfRays;
% %     end
% % end
% % 
% % wInitxRay(off:end) = 0;
% % wInitProt(1:sum(a(1:5))) = 0;
% % %wInitProt(sum(a(1:5))+1 : end ) = 0;
% % wInitCarb(:)=0;
% % doseXRay   = reshape(dij.physicalDose{1} * wInitxRay,[ct.cubeDim]);
% % doseProton = reshape(dij.physicalDose{1} * wInitProt,[ct.cubeDim]);
% % doseCarb = reshape(dij.physicalDose{1} * wInitCarb,[ct.cubeDim]);
% % doseTot    = reshape(dij.physicalDose{1} * wTot,[ct.cubeDim]);
% % 
% % slice = round(pln(1).isoCenter(1,3)/ct.resolution.z);
% % 
% % figure,
% % subplot(221),imagesc(doseXRay(:,:,slice)),colorbar
% % subplot(222),imagesc(doseProton(:,:,slice)),colorbar
% % subplot(223),imagesc(doseCarb(:,:,slice)),colorbar
% % subplot(224),imagesc(doseTot(:,:,slice)),colorbar
% 
% %%
% % figure, 
% % surf ( cube(:,:,slice))
% % %%
% % 
% % figure;
% % hold on,
% % plot(resultGUI.physicalDose(:,80,80));
% % plot(doseXRay(:,80,80));
% % plot(doseProton(:,80,80));
% % plot(doseCarb(:,80,80));
% % legend('total', 'Xray','Proton', 'Carbon');
% % title ('Default constraints, only Physical dose');
% % hold off;
% 
% 
% 
% %%

for i = 1:size(cst,1)
    if ~isfield(cst{i,5},'visibleColor')
        colorAssigned = false;
        break;
    elseif isempty(cst{i,5}.visibleColor)
        colorAssigned = false;
        break;
    end
end

% assign color if color assignment is not already present or inconsistent
if colorAssigned == false
  m         = 64;
  colorStep = ceil(m/size(cst,1));
  colors    = colorcube(colorStep*size(cst,1));
  % spread individual VOI colors in the colorcube color palette
  colors    = colors(1:colorStep:end,:);
  for i = 1:size(cst,1)
    cst{i,5}.visibleColor = colors(i,:);
  end
end

doseIsoLevels = [0.2:0.1:1.2]*2;
slice = round(pln(1).isoCenter(1,3)/ct.resolution.z);
c = jet(100);

% cube =  result_copy{1}.RBExD + result_copy{2}.RBExD;% + result_copy{3}.RBExD ;
% 
% 
% 
% % % 
% % % for i=(slice-20):(slice +20)
% % %     clf;
% % %     
% figure; 
% [hCMap,hDose,hCt,hContour,hIsoDose] =...
%     matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.RBExD,3,slice,[],[],[],[],[],doseIsoLevels);
% % % % pause(0.1);
% % % end
figure;

[hCMap,hDose,hCt,hContour,hIsoDose] =...
    matRad_plotSliceWrapper(gca,ct,cst,1,cube,3,slice,[],[],[],c,[],doseIsoLevels);
figure; 
[hCMap,hDose,hCt,hContour,hIsoDose] =...
    matRad_plotSliceWrapper(gca,ct,cst,1,result_copy{1}.RBExD,3,slice,[],[],[],c,[],doseIsoLevels);

figure; 
[hCMap,hDose,hCt,hContour,hIsoDose] =...
    matRad_plotSliceWrapper(gca,ct,cst,1,result_copy{2}.RBExD,3,slice,[],[],[],c,[],doseIsoLevels);
% 
% % 
