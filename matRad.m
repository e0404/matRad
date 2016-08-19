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

%  InputFolder = 'C:\MAtrad\data\4DCT\T6H_fuer_MB\biomech_samples';
%  numOfScen   = 2;
%  VOIs        = {'Blase','Haut','prostata_','Rektum','GTVPrimarius'};
%  [ct,cst]    = matRad_multScenImport(InputFolder,numOfScen,VOIs); 
 
load T6H_dose.mat
%load TKUH005_BPL.mat

%% multiple Scenarios
multScen.numOfCtScen         = ct.numOfCtScen; % number of imported ct scenarios
multScen.numOfShiftScen      = [0 0 0];        % number of shifts in x y and z direction       
multScen.shiftSize           = [3 3 3];     % equidistant: maximum shift [mm] / sampled: SD of normal distribution [mm]
multScen.shiftGenType        = 'equidistant';  % equidistant: equidistant shifts, sampled: sample shifts from normal distribution
multScen.shiftCombType       = 'individual';     % individual: no combination of shift scenarios, combined: combine shift scenarios
multScen.shiftGen1DIsotropy  = '+-';            % for equidistant shifts: '+-': positive and negative, '-': negative, '+': positive shift generation 
multScen.numOfRangeShiftScen = 0;              % number of absolute and/or relative range scnearios
%multScen.maxAbsRangeShift    = 0;              % maximum absolute over and undershoot in mm
%multScen.maxRelRangeShift    = 0;              % maximum relative over and undershoot in %
multScen.ScenCombType        = 'individual';   % individual: no combination of scenarios, allcombined: combine all scenarios
multScen                     = matRad_setMultScen(multScen);

%% meta information for treatment plan
pln.isoCenter       = matRad_getIsoCenter(cst,ct,0);
pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [0 90]; % [°]
pln.couchAngles     = [0 0]; % [Â°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.voxelDimensions = ct.cubeDim;
pln.radiationMode   = 'protons'; % either photons / protons / carbon
pln.bioOptimization = 'none'; % none: physical optimization; effect: effect-based optimization; RBExD: optimization of RBE-weighted dose
pln.numOfFractions  = 25;
pln.runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.machine         = 'HIT';
pln.minNrParticles  = 500000;
pln.LongitudialSpotSpacing = 3; %only relevant for HIT machine, not for generic
%% initial visualization and change objective function settings if desired
matRadGUI

%% generate steering file
stf = matRad_generateStf(ct,cst,pln,multScen);

%% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst,multScen);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst,multScen);
end

%% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% sequencing
if strcmp(pln.radiationMode,'photons') && (pln.runSequencing || pln.runDAO)
    %resultGUI = matRad_xiaLeafSequencing(resultGUI,stf,dij,5);
    resultGUI = matRad_engelLeafSequencing(resultGUI,stf,dij,5);
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

%% post processing
resultGUI = matRad_postprocessing(resultGUI, dij, pln, 25000000);

%% export Plan
matRad_export_HITXMLPlan_modified('T6H_stf_sv', 500000, 25000000, 'stfMode')  %500000 minNbParticles HIT Minimum für Patienten, minNrParticlesIES, scan path mode: 'stfMode', 'backforth','TSP' (very slow)

%% calc 4D dose
[resultGUI, delivery] = matRad_calc4dDose('T6H_stf_sv');  