%% Example Photon Treatment Plan with Direct aperture optimization
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
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
% (i) how to load patient data into matRad
% (ii) how to setup a photon dose calculation and 
% (iii) how to inversely optimize directly from command window in MatLab.
% (iv) how to apply a sequencing algorithm
% (v) how to run a direct aperture optimization
% (iv) how to visually and quantitatively evaluate the result

%% Patient Data Import
% Let's begin with a clear Matlab environment and import the head &
% neck patient into your workspace.
clc,clear,close all;
load('HEAD_AND_NECK.mat');

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.

pln.radiationMode   = 'photons';   % either photons / protons / carbon
pln.machine         = 'Generic';
pln.numOfFractions  = 30;

pln.propOpt.bioOptimization = 'none';    
pln.propStf.gantryAngles    = [0:72:359];
pln.propStf.couchAngles     = [0 0 0 0 0];
pln.propStf.bixelWidth      = 5;
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

%%
% Enable sequencing and direct aperture optimization (DAO).
pln.propOpt.runSequencing = 1;
pln.propOpt.runDAO        = 1;

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
% Lets generate dosimetric information by pre-computing dose influence 
% matrices for unit beamlet intensities. Having dose influences available 
% allows for subsequent inverse optimization.
dij = matRad_calcPhotonDose(ct,stf,pln,cst);

%% Inverse Planning for IMRT
% The goal of the fluence optimization is to find a set of beamlet weights 
% which yield the best possible dose distribution according to the 
% predefined clinical objectives and constraints underlying the radiation 
% treatment. Once the optimization has finished, trigger once the GUI to
% visualize the optimized dose cubes.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
matRadGUI;

%% Sequencing
% This is a multileaf collimator leaf sequencing algorithm that is used in 
% order to modulate the intensity of the beams with multiple static 
% segments, so that translates each intensity map into a set of deliverable 
% aperture shapes.
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,5);

%% DAO - Direct Aperture Optimization
% The Direct Aperture Optimization is an optimization approach where we 
% directly optimize aperture shapes and weights.
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);

%% Aperture visualization
% Use a matrad function to visualize the resulting aperture shapes
matRad_visApertureInfo(resultGUI.apertureInfo);

%% Indicator Calculation and display of DVH and QI
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI);
