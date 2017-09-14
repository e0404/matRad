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
load PROSTATE.mat
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
%pln.numOfFractions  = 38;
pln.runSequencing   = true; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = true; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.VMAT            = false; % 1/true: run VMAT, 0/false: don't
pln.dynamic         = true;
pln.halfFluOpt      = false; % indicates if you want to constrain half of each field to have 0 fluence (other half is compensated on the other side)
pln.machine         = 'Generic';


%% For VMAT
pln.runSequencing   = true;
pln.runDAO          = true;
pln.VMAT            = true;

pln.scaleDRx        = true;
pln.scaleDij        = true;
pln.jacobi          = true;
pln.dynamic         = true;


pln.numApertures = 7; %max val is pln.maxApertureAngleSpread/pln.minGantryAngleRes
pln.numLevels = 7;

pln.minGantryAngleRes = 4; %Bzdusek
pln.maxApertureAngleSpread = 28; %should be an even multiple of pln.minGantryAngleRes; Bzdusek

pln = matRad_VMATGantryAngles(pln,'new');


pln.gantryRotCst = [0 6]; %degrees per second
pln.defaultGantryRot = max(pln.gantryRotCst); %degrees per second
pln.leafSpeedCst = [0 6]*10; %mm per second
pln.defaultLeafSpeed = pln.leafSpeedCst(2);
pln.doseRateCst = [75 600]/60; %MU per second
pln.defaultDoseRate = pln.doseRateCst(2);
pln.maxLeafTravelPerDeg = pln.leafSpeedCst(2)/pln.defaultGantryRot;

pln.halfFluOpt = false;
pln.halfFluOptMargin = 10; % mm

recalc.doRecalc = 0;
recalc.dynamic = 0;
recalc.interpNew = 0;
recalc.pln = pln;
recalc.pln.minGantryAngleRes = 1;

%pln.halfFluOpt = false;



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
%dij.weightToMU = 100*(100/90)^2*(67/86)*(110/105)^2*(90/95)^2;
dij.weightToMU = 100;
dij.scaleFactor = 1;

%this is equal to multiplication of factors:
% - factor when reference conditions are equal to each other (100)
% - inverse square factor to get same SSD
% - PDD factor (evaluated at SSD = 100 cm) (Podgorsak IAEA pg. 183)
% - Mayneord factor to move SSD from 100 cm to 85 cm

%At TOH: 100 cm SAD, 5 cm depth, 10x10cm2
%At DKFZ: 95 cm SAD, 10 cm depth, 10x10cm2

%% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf,0);

%% sequencing
if strcmp(pln.radiationMode,'photons') && (pln.runSequencing || pln.runDAO)
    %resultGUI = matRad_xiaLeafSequencing(resultGUI,stf,dij,5);
    %resultGUI = matRad_engelLeafSequencing(resultGUI,stf,dij,5);
    resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);
    %resultGUI = matRad_svenssonLeafSequencing(resultGUI,stf,dij,pln,0);
end

%% DAO
if strcmp(pln.radiationMode,'photons') && pln.runDAO
   resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);
   %matRad_visApertureInfo(resultGUI.apertureInfo);
   
   recalc = matRad_doseRecalc(dij,cst,pln,recalc,ct,resultGUI.apertureInfo);
end

%% start gui for visualization of result
matRadGUI

%% dvh
matRad_calcDVH(resultGUI,cst,pln)

