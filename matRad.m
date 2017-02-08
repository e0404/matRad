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

clear
close all
clc

% load patient data, i.e. ct, voi, cst

%load HEAD_AND_NECK
%load TG119.mat
%load PROSTATE.mat
%load LIVER.mat
%load BOXPHANTOM.mat
%load([pwd filesep 'patients' filesep 'BOXPHANTOM_TINY.mat']);
%load(['/Volumes/WS_exFat/TG119/VWWC/TG119_VWWC.mat']);

load('/Volumes/WS_exFat/TG119/verification/TG119.mat')
%cst{2,5}.alphaX = 0.1; cst{2,5}.betaX = 0.01;
% 
% cst{2,6}(2,1)            = cst{2,6}(1);
% cst{2,6}(2,1).robustness = 'VWWC';

pln.exportInfluenceDataToASCII = false;

%% initial visualization and change objective function settings if desired
%matRadGUI

%% meta information for treatment plan
pln.isoCenter       = matRad_getIsoCenter(cst,ct,0);
pln.bixelWidth      = 4; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [0 45 315]; %[0:72:359]; % [°]
pln.couchAngles     = [0 0 0]; % [Â°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.voxelDimensions = ct.cubeDim;
pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.bioOptimization = 'MCN_RBExD';   % none: physical optimization;                                   const_RBExD; constant RBE of 1.1;  
                                     % LSM_effect;  variable RBE Linear Scaling Model (effect based); LSM_RBExD;  variable RBE Linear Scaling Model (RBExD based)
                                     % LEMIV_effect: effect-based optimization;                       LEMIV_RBExD: optimization of RBE-weighted dose
pln.numOfFractions         = 1;
pln.runSequencing          = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO                 = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.machine                = 'GenericLET';%GenericLET
pln.minNrParticles         = 500000;
pln.LongitudialSpotSpacing = 3;      % only relevant for HIT machine, not generic

%% initial visualization and change objective function settings if desired
%matRadGUI

%% retrieve model parameters
pln = matRad_getBioModel(pln);

%% set plan uncertainties for robust optimization
[cst,pln] = matRad_setPlanUncertainties(ct,cst,pln);

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);

%% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst,false);
    %dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst,false);
end

if pln.exportInfluenceDataToASCII
   addpath([pwd filesep 'IO' filesep 'MDACC']);
   [ flagSuccess ] = matRad_exportInfluenceDataToASCII(cst,stf,pln,dij,['/Volumes/WS_exFat/patient1']);
end


%% inverse planning for imrt
diary('OptimizationOutput')
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
diary off
%% sequencing
if strcmp(pln.radiationMode,'photons') && (pln.runSequencing || pln.runDAO)
    %resultGUI = matRad_xiaLeafSequencing(resultGUI,stf,dij,5);
    %resultGUI = matRad_engelLeafSequencing(resultGUI,stf,dij,5);
    resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,5);
end

%% DAO
if strcmp(pln.radiationMode,'photons') && pln.runDAO
   resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);
   matRad_visApertureInfo(resultGUI.apertureInfo);
end

%% start gui for visualization of result
% addpath([pwd filesep 'internal' filesep 'utilities']);


resultGUI = matRad_getBeamContributions(resultGUI,cst,stf,dij,'physicalDose');
resultGUIMDACC = matRad_getBeamContributions(resultGUIMDACC,cst,stf,dij,'physicalDose');
% 
% %matRadGUI
% 
% %% dvh
% matRad_calcDVH(resultGUI,cst,pln)

resultGUI = matRad_calcCubes(resultGUI.w,dij,cst,1);

resultGUIMDACC = matRad_calcCubes(VarName2,dij,cst,1);
