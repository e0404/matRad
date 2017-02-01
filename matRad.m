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
%load BOXPHANTOM_TINY.mat
%load(['/Volumes/WS_exFat/patient1/patient1_222mm.mat']);


%% initial visualization and change objective function settings if desired
%matRadGUI

pln.exportInfluenceDataToASCII = true;
%% meta information for treatment plan
pln.isoCenter       = matRad_getIsoCenter(cst,ct,0);
pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [0 45 315]; %[0:72:359]; % [°]
pln.couchAngles     = [0 0 0]; % [Â°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.voxelDimensions = ct.cubeDim;
pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.bioOptimization = 'LSM_RBExD';   % none: physical optimization;                                   const_RBExD; constant RBE of 1.1;  
                                     % LSM_effect;  variable RBE Linear Scaling Model (effect based); LSM_RBExD;  variable RBE Linear Scaling Model (RBExD based)
                                     % LEMIV_effect: effect-based optimization;                       LEMIV_RBExD: optimization of RBE-weighted dose
pln.robOpt          = 'none';        % none: optimize nominal scenario;     VWWC: voxel wise worst case optimization, 'probabilistic';  probabilistic optimization   
                                     % COWC: composite worst case;   coverage: coverage based optimization  
                                     % right now it doesnt allow different optimization techniques for different structures
pln.numOfFractions         = 25;
pln.runSequencing          = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO                 = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.machine                = 'GenericLET';
pln.minNrParticles         = 500000;
pln.LongitudialSpotSpacing = 3; %only relevant for HIT machine, not generic


%% initial visualization and change objective function settings if desired
%matRadGUI

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

%% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
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
matRadGUI

%% dvh
matRad_calcDVH(resultGUI,cst,pln)


