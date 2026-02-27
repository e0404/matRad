%% Example: Proton Treatment Plan with subsequent spot removal for Monte Carlo
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2021 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
%% set matRad runtime configuration
matRad_rc; %If this throws an error, run it from the parent directory first to set the paths

%% Patient Data Import
% Let's begin with a clear Matlab environment and import the prostate
% patient into your workspace

load('PROSTATE.mat');

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
pln.bioModel      = 'constRBE';
pln.multScen      = 'nomScen';

%%
% for particles it is possible to also calculate the LET disutribution
% alongside the physical dose. Therefore you need to activate the
% corresponding option during dose calculcation
pln.propDoseCalc.calcLET = 0;
                                       
%%
% Now we have to set the remaining plan parameters.
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 3;
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

%% Inverse Optimization for IMPT
% The goal of the fluence optimization is to find a set of bixel/spot 
% weights which yield the best possible dose distribution according to the 
% clinical objectives and constraints underlying the radiation treatment
pln.propOpt.quantityOpt = quantityOpt;
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Spot removal
% instantiate spot removal class
spotRemover = matRad_SpotRemovalDij(dij,resultGUI.w);

spotRemover.removalMode = 'relative';
spotRemover.propSpotRemoval.relativeThreshold = 0.05;
resultGUI2 = spotRemover.reoptimize(cst,pln);

% numOfRemovedSpots = sr_cfg.numOfRemovedSpots;

% stf2 = sr_cfg.getStf(stf);
% dij2 = sr_cfg.getDij;
% weight2 = sr_cfg.getWeights;
% weightLogical = sr_cfg.getLogical;

%% Plot difference of the doses
matRad_compareDose(resultGUI.RBExDose,resultGUI2.RBExDose,ct,cst);


