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
% Let's begin with a clear Matlab environment and import the TG119 patient
% into your workspace
matRad_rc

load TG119.mat

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.

% meta information for treatment plan
pln.numOfFractions  = 30;
pln.radiationMode   = 'photons';            % either photons / protons / helium / carbon / brachy
pln.machine         = 'Generic';            % generic for RT / LDR or HDR for BT

pln.bioModel = 'none';      % none: for photons, protons, carbon, brachy    % constRBE: constant RBE for photons and protons 
                            % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                            % LEM: Local Effect Model for carbon ions       % HEL: data-driven RBE parametrization for helium

pln.multScen = 'nomScen';   % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'         

% beam geometry settings
pln.propStf.bixelWidth      = 5;            % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.maxGantryAngleSpacing    = 30;   % [°] / max gantry angle spacing for dose calculation
pln.propStf.maxDAOGantryAngleSpacing = 60;  % [°] / max gantry angle spacing for DAO
pln.propStf.maxFMOGantryAngleSpacing = 180;  % [°] / max gantry angle spacing for FMO
pln.propStf.startingAngle = -180;           % [°] / starting angle for VMAT
pln.propStf.finishingAngle = 180;           % [°] / finishing angle for VMAT
pln.propStf.couchAngle      = 0;            % [°]
pln.propStf.isoCenter       = matRad_getIsoCenter(cst,ct,0);
pln.propStf.generator       = 'PhotonVMAT';
pln.propStf.continuousAperture  = false;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

% sequencing settings
pln.propSeq.runSequencing   = true;  % true: run sequencing, false: don't / will be ignored for particles and also triggered by runDAO below
pln.propSeq.sequencer       = 'siochi';
pln.propSeq.numLevels       = 7;

% optimization settings
pln.propOpt.quantityOpt         = 'physicalDose';   % Quantity to optimizer (could also be RBExDose, BED, effect)
pln.propOpt.optimizer           = 'IPOPT';          % We can also utilize 'fmincon' from Matlab's optimization toolbox
pln.propOpt.runDAO              = true;             % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runVMAT             = true;
pln.propOpt.preconditioner      = true;

%pln.propOpt.VMAToptions.machineConstraintFile = [pln.radiationMode '_' pln.machine];

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
% Lets generate dosimetric information by pre-computing dose influence 
% matrices for unit beamlet intensities. Having dose influences available 
% allows for subsequent inverse optimization.
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

%% Inverse Planning for IMRT
% The goal of the fluence optimization is to find a set of beamlet weights 
% which yield the best possible dose distribution according to the 
% predefined clinical objectives and constraints underlying the radiation 
% treatment. In VMAT, FMO is done only at the angles in the
% FMOGantryAngles set. Once the optimization has finished, trigger once the GUI to
% visualize the optimized dose cubes.
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);
matRadGUI;

%% Sequencing
% This is a multileaf collimator leaf sequencing algorithm that is used in 
% order to modulate the intensity of the beams with multiple static 
% segments, so that translates each intensity map into a set of deliverable 
% aperture shapes. The fluence map at each angle in the initGantryAngles
% set is sequenced, with the resulting apertures spread to neighbouring
% angles from the optGantryAngles set.
resultGUI = matRad_sequencing(resultGUI,stf,dij,pln);

%% DAO - Direct Aperture Optimization
% The Direct Aperture Optimization is an optimization approach where we 
% directly optimize aperture shapes and weights at the angles in the
% optGantryAngles set.  The gantry angle speed, leaf speed, and MU rate are
% constrained by the min and max values specified by the user.
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);

%% Aperture visualization
% Use a matrad function to visualize the resulting aperture shapes
matRad_visApertureInfo(resultGUI.apertureInfo);

%% Indicator Calculation and display of DVH and QI
resultGUI = matRad_planAnalysis(resultGUI,ct,cst,stf,pln);

%% Calculate delivery metrics

resultGUI = matRad_calcDeliveryMetrics(resultGUI,pln,stf);

