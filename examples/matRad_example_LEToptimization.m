%% Example: Proton LET optimization
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
% for particles with tabulated LET it is possible to also calculate the LET 
% disutribution alongside the physical dose. Therefore you need to activate
% the corresponding option during dose calculcation
pln.propDoseCalc.calcLET = true;
                                       
%%
% Now we have to set the remaining plan parameters.
pln.numOfFractions        = 1;
pln.propStf.gantryAngles  = [0 120 240];
pln.propStf.couchAngles   = [0 0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.longitudinalSpotSpacing = 3;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
% Lets generate dosimetric information by pre-computing dose influence 
% matrices for unit beamlet intensities. Having dose influences available 
% allows for subsequent inverse optimization. 
dij = matRad_calcParticleDose(ct,stf,pln,cst);


%% Run First Fluence Optimization without LET objectives
% The goal of the fluence optimization is to find a set of bixel/spot 
% weights which yield the best possible dose distribution according to the 
% clinical objectives and constraints underlying the radiation treatment
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot dose and LET
plotF = figure;
pDose = subplot(2,2,1);
matRad_plotSliceWrapper(pDose,ct,cst,1,resultGUI.RBExDose,3,65);

pLET = subplot(2,2,2);
matRad_plotSliceWrapper(pLET,ct,cst,1,resultGUI.LET,3,65);


%% Add LET objectives
%cst{1,6}{2} = struct(LETObjectives.matRad_MeanLET(1000));
cst{1,6}{2} = struct(LETObjectives.matRad_SquaredOverLET(100,0));

cst{2,6}{2} = struct(LETObjectives.matRad_SquaredUnderLET(100,50));
%cst{2,6}{2} = struct(LETObjectives.matRad_SquaredOverLET(100,0));
%cst{2,6}{2} = struct(LETObjectives.matRad_SquaredDeviationLET(100,0));

%% Second Optimization with LET objectives
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot dose and LET
%plotF = figure;
pDose = subplot(2,2,3);
matRad_plotSliceWrapper(pDose,ct,cst,1,resultGUI.RBExDose,3,65);

pLET = subplot(2,2,4);
matRad_plotSliceWrapper(pLET,ct,cst,1,resultGUI.LET,3,65);

%% Inverse Optimization for IMPT
matRadGUI




