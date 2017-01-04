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
load TG119.mat
%load PROSTATE.mat
%load LIVER.mat
%load BOXPHANTOM.mat

% meta information for treatment plan
pln.isoCenter       = matRad_getIsoCenter(cst,ct,0);
pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [0:72:359]; % [°]
pln.couchAngles     = [0 0 0 0 0]; % [°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.voxelDimensions = ct.cubeDim;
pln.radiationMode   = 'photons';     % either photons / protons / carbon
pln.bioOptimization = 'none';        % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                     % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.numOfFractions  = 30;
pln.runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.VMAT            = false; % 1/true: run VMAT, 0/false: don't
pln.machine         = 'Generic';


%% For VMAT
pln.runSequencing   = true;
pln.runDAO          = true;
pln.VMAT            = true;

pln.numApertures = 7;
pln.numLevels = 3;

pln.minGantryAngleRes = 4;
pln.maxApertureAngleSpread = 20; %should be 1/2
%Why should this be smaller than 10\deg?
pln.numInitGantryAngles = max([360/pln.maxApertureAngleSpread 360/(pln.numApertures*pln.minGantryAngleRes)]);

pln.initGantryAngleSpacing = 360/pln.numInitGantryAngles;
pln.initGantryAngles = pln.initGantryAngleSpacing/2+pln.initGantryAngleSpacing*(0:(pln.numInitGantryAngles-1)); %pln.optGantryAngleSpacing*(pln.numApertures-1)/2+pln.initGantryAngleSpacing*(0:(pln.numInitGantryAngles-1));

pln.optGantryAngleSpacing = pln.initGantryAngleSpacing/pln.numApertures;
pln.optGantryAngles = pln.optGantryAngleSpacing/2+pln.optGantryAngleSpacing*(0:(pln.numInitGantryAngles*pln.numApertures-1)); %0:pln.optGantryAngleSpacing:360;

pln.optToGantryAngleSpacingFactor = ceil(pln.optGantryAngleSpacing/4);

pln.gantryAngleSpacing = pln.optGantryAngleSpacing/pln.optToGantryAngleSpacingFactor; %ideally should be spaced every 2 or 4 degrees; gantry spacing that dij is performed
pln.gantryAngles    = pln.gantryAngleSpacing/2+pln.gantryAngleSpacing*(0:(pln.numInitGantryAngles*pln.numApertures*pln.optToGantryAngleSpacingFactor-1)); %0:pln.gantryAngleSpacing:360;
pln.couchAngles     = 0*pln.gantryAngles;

pln.numOfBeams      = numel(pln.gantryAngles);

pln.gantryRotCst = [0 6]; %degrees per second
pln.defaultGantryRot = mean(pln.gantryRotCst); %degrees per second
pln.leafSpeedCst = [-inf inf]; %cm per second
pln.doseRateCst = [0 inf]; %MU per second


%% initial visualization and change objective function settings if desired
matRadGUI

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);

%% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
end
dij.weightToMU = 100;
%100 cm SAD, 5 cm depth, 10x10cm2

%% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% sequencing
if strcmp(pln.radiationMode,'photons') && (pln.runSequencing || pln.runDAO)
    %resultGUI = matRad_xiaLeafSequencing(resultGUI,stf,dij,5);
    %resultGUI = matRad_engelLeafSequencing(resultGUI,stf,dij,5);
    resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln.numLevels,0,pln.VMAT,pln.numApertures);
    %resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,7,0,0,0);
end

%% DAO
if strcmp(pln.radiationMode,'photons') && pln.runDAO
   resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);
   matRad_visApertureInfo(resultGUI.apertureInfo);
end

%% start gui for visualization of result
matRadGUI

%% dvh
matRad_calcDVH(resultGUI,cst,pln)

