%% 4D dose calculation workflow
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In this example we will show 
% (i) the structure of 4D data within matRad
% (ii) how to perform standard treatment planning
% (iii) how to run a dose recalculation considering interplay effects 
 
%% set matRad runtime configuration
matRad_rc

%% Load data, add generic 4D information, and display 'moving' geometry
load BOXPHANTOM.mat

%%

amplitude    = [0 3 0]; % [voxels]
numOfCtScen  = 5;
motionPeriod = 2.5; % [s] 

[ct,cst] = matRad_addMovement(ct, cst,motionPeriod, numOfCtScen, amplitude);
% Set up a plan, compute dose influence on all phases, conventional optimization
% meta information for treatment plan
pln.numOfFractions  = 30;
pln.radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln.machine         = 'Generic';

% beam geometry settings
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longitudinalSpotSpacing = 5;      % only relevant for HIT machine, not generic
pln.propStf.gantryAngles    = [90]; 
pln.propStf.couchAngles     = [0]; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

%optimization settings
pln.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

quantityOpt  = 'RBExD';     % options: physicalDose, effect, RBExD
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions

scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen' 

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,scenGenType);

%%
% generate steering file
stf = matRad_generateStf(ct,cst,pln);

%% 
% dose calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

%% 
% inverse planning for imrt on a static CT
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% 
% post processing
% This step is necessary to remove beam spots with too few particles that
% cannot not be delivered, dose is recalculated accordingly
resultGUI = matRad_postprocessing(resultGUI, dij, pln, cst, stf) ; 

%%
% calc 4D dose
% make sure that the correct pln, dij and stf are loeaded in the workspace
[resultGUI, timeSequence] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI); 

% plot the result in comparison to the static dose
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z); 

figure 

subplot(2,2,1)
imagesc(resultGUI.RBExD(:,:,slice)),colorbar, colormap(jet);
title('static dose distribution [Gy (RBE)]')
axis equal

subplot(2,2,3)
imagesc(resultGUI.accRBExD(:,:,slice)),colorbar, colormap(jet); 
title('accumulated (4D) dose distribution [Gy (RBE)]')
axis equal

subplot(2,2,2)
imagesc(resultGUI.RBExD(:,:,slice) - resultGUI.accRBExD(:,:,slice)) ,colorbar, colormap(jet); 
title('static dose distribution - accumulated (4D) dose distribution [Gy (RBE)]')

axis equal

