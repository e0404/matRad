%% MatRad DAO example
%
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

%% Inizialising
% After loading the patient ct and cst, we proceed setting the plan of the 
% treatment. 
% We define, in this example, a 9-rays-treatment with photons and enable the
% option 'pln.runDAO'.
clear
close all
clc

% load patient data, i.e. ct, voi, cst
load PROSTATE.mat

% meta information for treatment plan
pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [0:40:359]; % [°]
pln.couchAngles     = [0 0 0 0 0 0 0 0 0]; % [°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.isoCenter       = ones(pln.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.voxelDimensions = ct.cubeDim;
pln.radiationMode   = 'photons';     % either photons / protons / carbon
pln.bioOptimization = 'none';        % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                     % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.numOfFractions  = 30;
pln.runSequencing   = true; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = true; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.machine         = 'Generic';

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);

%% dose calculation
dij = matRad_calcPhotonDose(ct,stf,pln,cst);

%% inverse planning for imrt
% This function optimizes the fluence of the beam, giving back the weights
% used to scale the number of photons in every ray of the beam, improving
% the accuracy of the simulation
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% sequencing
% This is a multileaf collimator leaf sequencing algorithm that is used in 
% order to modulate the intensity of the beams with multiple static 
% segments, so that translates each intensity map into a set of deliverable 
% aperture shapes; according to Siochi (1999).
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,5);

%% DAO
% The Direct Aperture Optimization is an automated planning system, only
% possibble for photons in which we bypass the traditional intensity
% optimization, and instead directly optimize the shapes and the weights of
% the apertures. This technique allows the user to specify the maximum
% number of apertures per beam direction, and hence provides significant
% control over the complexity of the treatment delivery.
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);
matRad_visApertureInfo(resultGUI.apertureInfo);

%% indicator calculation
% Here we call the indicator wrapper function that is just a wrapper
% function to calculate multiple quality indicators like D98, D95. One
% basically calculates default quality indicators for a given dose distribution.
% If you don't see clearly the images, go in the folder 'html' in your
% working folder and you should look at the images, starting with the same
% name of this file, numbered from 3 to 11.
cst = matRad_indicatorWrapper(cst,pln,resultGUI);

%% show DVH and QI
% This last function allow you to plot the Dose Volume Histogram and the
% dose statistics of all the VOIs in side the ct.
matRad_showDVH(cst,pln)

