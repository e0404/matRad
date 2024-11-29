%% 4D dose calculation workflow
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
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

amplitude    = [0 3 0]; % [voxels]
numOfCtScen  = 5;
motionPeriod = 2.5; % [s] 

[ct,cst] = matRad_addMovement(ct, cst,motionPeriod, numOfCtScen, amplitude,'dvfType','pull');

% Set up a plan, compute dose influence on all phases, conventional optimization
% meta information for treatment plan
pln.numOfFractions  = 30;
pln.radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln.machine         = 'Generic';
pln.bioModel        = 'constRBE';

% retrieve scenarios for dose calculation and optimziation.
% A nominal Scenario Model will consider all 4D scenarios as they are not
% "uncertainty" per se
pln.multScen = matRad_multScen(ct, 'nomScen');

% beam geometry settings
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longitudinalSpotSpacing = 5;      % only relevant for HIT machine, not generic
pln.propStf.gantryAngles    = [90]; 
pln.propStf.couchAngles     = [0]; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);


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
slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
slice = slice(3);

figure 

subplot(2,2,1)
imagesc(resultGUI.RBExDose(:,:,slice)),colorbar, colormap(jet);
title('static dose distribution [Gy (RBE)]')
axis equal

subplot(2,2,3)
imagesc(resultGUI.accRBExDose(:,:,slice)),colorbar, colormap(jet); 
title('accumulated (4D) dose distribution [Gy (RBE)]')
axis equal

subplot(2,2,2)
imagesc(resultGUI.RBExDose(:,:,slice) - resultGUI.accRBExDose(:,:,slice)) ,colorbar, colormap(jet); 
title('static dose distribution - accumulated (4D) dose distribution [Gy (RBE)]')

axis equal

