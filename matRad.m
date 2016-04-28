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

% InputFolder = 'E:\Mescher\13_BIOM_model\01_BioMechModel\Input';
% numOfScen   = 3;
% VOIs        = {'Blase','Haut','prostata_','Rektum','GTVPrimarius'};
% [ct,cst]    = matRad_multScenImport(InputFolder,numOfScen,VOIs); 
% load T6H.mat

%% multiple Scenarios
multScen.numOfCtScen         = ct.numOfCtScen; % number of imported ct scenarios
multScen.numOfShiftScen      = [9 9 9];        % number of shifts in x y and z direction       
multScen.shiftSize           = [9 9 9];        % equidistant: maximum shift [mm] / sampled: SD of normal distribution [mm]
multScen.shiftGenType        = 'equidistant';  % equidistant: equidistant shifts, sampled: sample shifts from normal distribution
multScen.shiftGen1DIsotropy  = '+';            % for equidistant shifts: '+-': positive and negative, '-': negative, '+': positive shift generation 
multScen.shiftCombType       = 'combined';     % for equidistant shifts: combined, individual
multScen.numOfRangeShiftScen = 0;              % number of absolute and/or relative range scnearios
multScen.maxAbsRangeShift    = 0;              % maximum absolute over and undershoot in mm
multScen.maxRelRangeShift    = 0;              % maximum relative over and undershoot in %
multScen.ScenCombType        = 'individual';   % individual: no combination of scenarios, allcombined: combine all scenarios
multScen                     = matRad_setMultScen(multScen);

%% coverage based cst manipulation
%load('E:\Mescher\15_DCH_objectiv_tests\01_PROSTATE_InisotropicShifts\02_PROSTATE_5photonBeams_10XYZshifts_ProbRing_5mmBixel\standard_prob_CST')
load('E:\Mescher\15_DCH_objectiv_tests\01_PROSTATE_InisotropicShifts\02_PROSTATE_5photonBeams_10XYZshifts_ProbRing_5mmBixel\standard_nonprob_CST')

cst = matRad_coverageBasedCstManipulation(cst,ct,multScen,'probWeighting');

%% meta information for treatment plan
pln.isoCenter       = matRad_getIsoCenter(cst,ct,0);
pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [0:72:359]; % [°]
pln.couchAngles     = [0 0 0 0 0]; % [Â°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.voxelDimensions = ct.cubeDim;
pln.radiationMode   = 'photons'; % either photons / protons / carbon
pln.bioOptimization = 'none'; % none: physical optimization; effect: effect-based optimization; RBExD: optimization of RBE-weighted dose
pln.numOfFractions  = 1;
pln.runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.machine         = 'Generic';

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
resultGUI = matRad_fluenceOptimization(dij,ct,cst,pln,multScen);

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

