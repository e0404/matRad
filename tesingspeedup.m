%% NEXT STEP IS TO ITERATE AND LOG TO FIND THE BEST TRADE OFF VALUES
%  ALSO TO MAKE A PLOT :P
% iterate over sampling percentages
% save : ogt, spt, gamma , voxel reduction

% ogt: original timing
% spt: sparcesampled timing

% load patient data, i.e. ct, voi, cst
clc
clear all

matRad_rc

 load TG119.mat
% load PROSTATE.mat
% load LIVER.mat
% load BOXPHANTOM.mat


%% meta information for treatment plan (1)
pln(1).numOfFractions  = 10;
pln(1).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln(1).machine         = 'Generic';

% beam geometry settings
pln(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(1).propStf.gantryAngles    = [0:72:359]; % [?] ;
pln(1).propStf.couchAngles     = zeros(numel(pln(1).propStf.gantryAngles),1); % [?] ;
pln(1).propStf.numOfBeams      = numel(pln(1).propStf.gantryAngles);
pln(1).propStf.isoCenter       = ones(pln(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(1).propOpt.spatioTemp      = 0;
pln(1).propOpt.STScenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln(1).propDoseCalc.doseGrid.resolution = ct.resolution; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution.y = 3; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution.z = 3; % [mm]

quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons
% MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons
% LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'

% retrieve bio model parameters
pln(1).bioParam = matRad_bioModel(pln(1).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln(1).multScen = matRad_multScen(ct,scenGenType);
%% generate steering file
stf = matRad_generateStf(ct,cst,pln);

%% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij1 = matRad_calcPhotonDose(ct,stf,pln,cst,param);
    %dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'helium') || strcmp(pln.radiationMode,'carbon')
    dij1 = matRad_calcParticleDose(ct,stf,pln,cst,param);
end

%% Resampled Dij
% this could be moved into matRad_fluencOptimization too

margin_size = 2;            % in number of voxels
sparCity = [ 0.7 0.1];        % fraction of voxels sampled in [target OAR]
[cst,dij2] = matRad_dijSampler(cst,ct,dij1,margin_size, sparCity);

% number of voxels in first Dij
ognum = numel(find(sum(dij1.physicalDose{1},2)));
% number of voxels in second Dij
spnum = numel(find(sum(dij2.physicalDose{1},2)));
totalfxnVox = spnum/ognum 

%%  running the OG one
tic
result1  = matRad_fluenceOptimization(dij1,cst,pln,param);
ogt = toc;
cst1 = cst;
%% running faster one
tic
res_hold =  matRad_fluenceOptimization(dij2,cst,pln,param);
result2 = matRad_calcCubes(res_hold.w,dij1);  % dose calculation to always be done on full matrix
spt = toc;
cst = cst1;
% matRad_gammaIndex(cube1,cube2,resolution,criteria,slice,n,localglobal,cst)
gammaCube = matRad_gammaIndex(result1.physicalDose,result2.physicalDose,[3 3 3],[1 1],80);

% gamma = matRad_gammaIndex(
% [gammaCube,gammaPassRate,hfig] = matRad_compareDose(result1.physicalDose,result2.physicalDose, ct, cst);
% gammaPassRate