%% Example: Proton XBDLET optimization
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% In this example we will show 


%% Patient Data Import
% Let's begin with a clear Matlab environment and import the prostate
% patient into your workspace

matRad_rc; %If this throws an error, run it from the parent directory first to set the paths
clear all;
load('TG119.mat');

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines the most 
% important cornerstones of your treatment plan.

%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon. In this
% example we would like to use protons for treatment planning. Next, we
% need to define a treatment machine to correctly load the corresponding 
% base data. matRad features generic base data in the file
% 'proton_Generic.mat'; consequently the machine has to be set accordingly
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

%%
% Define the flavor of biological optimization for treatment planning along
% with the quantity that should be used for optimization. Possible values 
% are (none: physical optimization; const_RBExD: constant RBE of 1.1; 
% LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of 
% RBE-weighted dose. As we use protons, we follow here the clinical 
% standard and use a constant relative biological effectiveness of 1.1. 
% Therefore we set bioOptimization to const_RBExD
pln.propOpt.bioOptimization = 'const_RBExD';

                                       
%%
% Now we have to set the remaining plan parameters.
pln.numOfFractions        = 1;
pln.propStf.gantryAngles  = [0];
pln.propStf.couchAngles   = [0];
pln.propStf.bixelWidth    = 5;
pln.propStf.longitudinalSpotSpacing = 3;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

pln.propDoseCalc.calcLET = true;

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
% Lets generate dosimetric information by pre-computing dose influence 
% matrices for unit beamlet intensities. Having dose influences available 
% allows for subsequent inverse optimization.
dij = matRad_calcParticleDose(ct,stf,pln,cst);

%% Add the factor for XBD
dij.c = 0.04;

%% Add the dose rate / MU settings
%Set a fixed current in the dij
dij.fixedCurrent = 300;

%Set some MU data
dij.minMU = 300; % minimum MU to be produced
dij.MU = 5; %factor MU -> w (now means 1 MU = 0.16 * w = 0.16*1e6 particles)

%% Prescriptions for single fraction treatment
cst{1,6}{1}.parameters{1} = 2;
cst{2,6}{1}.parameters{1} = 5;
cst{3,6}{1}.parameters{1} = 2.5;

%% Run First Fluence Optimization without LET objectives
% The goal of the fluence optimization is to find a set of bixel/spot 
% weights which yield the best possible dose distribution according to the 
% clinical objectives and constraints underlying the radiation treatment
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot dose and LET
plotF = figure;
pDose = subplot(2,3,1);
matRad_plotSliceWrapper(pDose,ct,cst,1,resultGUI.RBExDose,3,65);
title(pDose,'RBExD [1.1]');

pDose = subplot(2,3,2);
matRad_plotSliceWrapper(pDose,ct,cst,1,resultGUI.XBD_LET,3,65);
title(pDose,'XBD_LET');

pLET = subplot(2,3,3);
matRad_plotSliceWrapper(pLET,ct,cst,1,resultGUI.LET,3,65);
title(pDose,'LET');

%%
% Now compare with XBD_LET optimization
pln.propOpt.bioOptimization = 'XBD_LET';
FLASH_doseRate = 40; % Gy/s
cst{1,6}{2} = struct(DADRObjectives.matRad_SquaredUnderDADR(1e6, FLASH_doseRate));

resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot
figure(plotF);
pDose = subplot(2,3,4);
matRad_plotSliceWrapper(pDose,ct,cst,1,resultGUI.RBExDose,3,65);
title(pDose,'RBExD [1.1]');

pDose = subplot(2,3,5);
matRad_plotSliceWrapper(pDose,ct,cst,1,resultGUI.XBD_LET,3,65);
title(pDose,'XBD_LET');

pLET = subplot(2,3,6);
matRad_plotSliceWrapper(pLET,ct,cst,1,resultGUI.LET,3,65);
title(pDose,'LET');

%% Inverse Optimization for IMPT
matRadGUI




