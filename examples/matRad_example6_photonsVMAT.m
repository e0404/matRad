%% Example Photon Treatment Plan with VMAT direct aperture optimization
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
% (ii) how to input necessary parameters in the pln structure
% (iii) how to setup a photon dose calculation
% (iv) how to inversely optimize fluence directly from command window in MatLab.
% (v) how to apply a sequencing algorithm
% (vi) how to run a VMAT direct aperture optimization
% (vii) how to visually and quantitatively evaluate the result

%% Patient Data Import
% Let's begin with a clear Matlab environment and import the head &
% neck patient into your workspace.
matRad_rc

load('HEAD_AND_NECK.mat');

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.

pln.radiationMode   = 'photons';   % either photons / protons / carbon
pln.machine         = 'Generic';

pln.numOfFractions  = 30;

% beam geometry settings
pln.propStf.bixelWidth = 10;

% optimization settings
pln.propOpt.bioOptimization = 'none';
pln.propOpt.runVMAT         = true;
pln.propOpt.runDAO          = true;
pln.propOpt.runSequencing   = true;
pln.propOpt.preconditioner  = true;
pln.propOpt.numLevels       = 7;
 
pln.propOpt.VMAToptions.machineConstraintFile = [pln.radiationMode '_' pln.machine];

pln.propOpt.VMAToptions.maxGantryAngleSpacing    = 2;      % Max gantry angle spacing for dose calculation
pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing = 4;      % Max gantry angle spacing for DAO
pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing = 28;     % Max gantry angle spacing for FMO

pln.propOpt.VMAToptions.startingAngle = 0;
pln.propOpt.VMAToptions.finishingAngle = 359;
pln.propOpt.VMAToptions.continuousAperture = 0;

pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

%%
% Generate dose calculation, DAO, and FMO angles from the parameters input
% above. FMO is performed only on the initGantryAngles set. In the DAO
% step, weights and leaf positions are optimized at the angles in the
% optGantryAngles set. Weights and leaf positions are interpolated at the
% angles in the gantryAngles set to increase the accuracy of the dose
% calculation (each iteration).

% FMO: optimize fluence on coarse subset of gantry angles
% Sequencing: select subset of apertures and spread to finer angles
% DAO: constrain for leaf speed, gantry rotation speed and MU rate

pln = matRad_VMATGantryAngles(pln,cst,ct);

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
% treatment. In VMAT, FMO is done only at the angles in the
% initGantryAngles set. Once the optimization has finished, trigger once the GUI to
% visualize the optimized dose cubes.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
matRadGUI;

%% Sequencing
% This is a multileaf collimator leaf sequencing algorithm that is used in 
% order to modulate the intensity of the beams with multiple static 
% segments, so that translates each intensity map into a set of deliverable 
% aperture shapes. The fluence map at each angle in the initGantryAngles
% set is sequenced, with the resulting apertures spread to neighbouring
% angles from the optGantryAngles set.
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);

%% DAO - Direct Aperture Optimization
% The Direct Aperture Optimization is an optimization approach where we 
% directly optimize aperture shapes and weights at the angles in the
% optGantryAngles set.  The gantry angle speed, leaf speed, and MU rate are
% constrained by the min and max values specified by the user.
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);

%% Aperture visualization
% Use a matrad function to visualize the resulting aperture shapes
matRad_visApertureInfo(resultGUI.apertureInfo);

%% Indicator Calculation and display of DVH and QI
cst = matRad_indicatorWrapper(cst,pln,resultGUI);
matRad_showDVH(cst,pln);

