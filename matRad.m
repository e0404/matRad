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
load LIVER.mat
%load BOXPHANTOM.mat


%% multiple Scenarios
multScen.numOfCtScen          = ct.numOfCtScen; % number of imported ct scenarios

multScen.numOfIntSegShiftScen = 0; %1000;              % number of internal segmentation shift scnearios     

multScen.numOfShiftScen       = [2 2 0];        % number of shifts in x y and z direction       
multScen.shiftSize            = [5 5 5];        % maximum shift [mm]
multScen.shiftSD              = [3 3 3];        % SD of normal distribution [mm]
multScen.shiftGenType         = 'equidistant';  % equidistant: equidistant shifts, sampled: sample shifts from normal distribution
multScen.shiftCombType        = 'individual';   % individual: no combination of shift scenarios, combined: combine shift scenarios, allcombined: create every possible shift combination
multScen.shiftGen1DIsotropy   = '+-';           % for equidistant shifts: '+-': positive and negative, '-': negative, '+': positive shift generation 

multScen.numOfRangeShiftScen  = 2;              % number of absolute and/or relative range scnearios
multScen.maxAbsRangeShift     = 5;              % maximum absolute over and undershoot in mm
multScen.maxRelRangeShift     = 0;              % maximum relative over and undershoot in %
multScen.ScenCombType         = 'individual';   % individual: no combination of scenarios, allcombined: combine all scenarios
multScen                      = matRad_setMultScen(multScen);

%% initial visualization and change objective function settings if desired
matRadGUI

%% coverage based cst manipulation
cst = matRad_coverageBasedCstManipulation(cst,ct,multScen,0,0);

%% meta information for treatment plan
pln.isoCenter       = matRad_getIsoCenter(cst,ct,0);
pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [0 90]; %[0:72:359]; % [°]
pln.couchAngles     = [0 0 ]; % [Â°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.voxelDimensions = ct.cubeDim;
pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.bioOptimization = 'none';        % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                     % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.numOfFractions  = 30;
pln.runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.machine         = 'HIT';
pln.minNrParticles  = 500000;
pln.LongitudialSpotSpacing = 3; %only relevant for HIT machine, not generic
%% initial visualization and change objective function settings if desired
matRadGUI

%% generate steering file
stf = matRad_generateStf(ct,cst,pln,multScen);

%% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst,multScen,false);
    %dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst,multScen,false);
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


