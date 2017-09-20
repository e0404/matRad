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
load BOXPHANTOM.mat

% meta information for treatment plan
pln.bixelWidth      = 10; % [mm]
pln.gantryAngles    = 0; % [°]
pln.couchAngles     = 0; % [°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.isoCenter       = ones(pln.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.voxelDimensions = ct.cubeDim;
pln.radiationMode   = 'photons';     % either photons / protons / carbon
pln.bioOptimization = 'none';        % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                     % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.numOfFractions  = 30;
pln.runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.machine         = 'Generic';

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);

%% dose calculation
dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst,500,1);

%% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

