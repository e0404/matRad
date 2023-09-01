%% Example: Photon Treatment Plan using VMC++ dose calculation
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

%% In this example we will show 
% (i) how to load patient data into matRad
% (ii) how to setup a photon dose calculation based on the VMC++ Monte Carlo algorithm 
% (iii) how to inversely optimize the beamlet intensities directly from command window in MATLAB. 
% (iv) how to visualize the result

%% set matRad runtime configuration
matRad_rc %If this throws an error, run it from the parent directory first to set the paths

%% Patient Data Import
% Let's begin with a clear Matlab environment and import the boxphantom
% into your workspace. 
load('BOXPHANTOM.mat');

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines the most
% important cornerstones of your treatment plan.

pln.radiationMode           = 'photons';  
pln.machine                 = 'Generic';
pln.numOfFractions          = 30;
pln.propStf.gantryAngles    = [0];
pln.propStf.couchAngles     = [0];
pln.propStf.bixelWidth      = 10;
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% Enable sequencing and direct aperture optimization (DAO).
pln.propOpt.runSequencing   = 1;
pln.propOpt.runDAO          = 1;

quantityOpt    = 'physicalDose';                                     
modelName      = 'none';  

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen');
% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

%%
pln.propMC.engine = 'TOPAS';
pln.propMC.beamProfile = 'phasespace';
pln.propMC.externalCalculation = true;
pln.propMC.infilenames.phaseSpaceSourcePhotons = 'SIEMENS_PRIMUS_6mv_15x15';
%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
dij = matRad_calcPhotonDose(ct,stf,pln,cst);

%% Inverse Optimization for IMRT
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
%% Sequencing
% This is a multileaf collimator leaf sequencing algorithm that is used in 
% order to modulate the intensity of the beams with multiple static 
% segments, so that translates each intensity map into a set of deliverable 
% aperture shapes.
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,5,1);
[pln,stf] = matRad_aperture2collimation(pln,stf,resultGUI.sequencing,resultGUI.apertureInfo);
%% Aperture visualization
% Use a matrad function to visualize the resulting aperture shapes
matRad_visApertureInfo(resultGUI.apertureInfo)
%% Plot the Resulting Dose Slice
% Just let's plot the transversal iso-center dose slice
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
figure,
imagesc(resultGUI.physicalDose(:,:,slice)),colorbar, colormap(jet)

%% Dose Calculation
%stf.ray.energy = [6,6,6];
%stf.ray.weight = [stf.ray.shapes.weight];
pln.propMC.numHistories = 1e8;
resultGUI_MC = matRad_calcPhotonDoseMC(ct,stf,pln,cst,1);

%% readout
foldername = 'E:\Code\matRad\topas\MCrun\photons_Generic_01-09-23';
pln = matRad_cfg.getDefaultClass(pln,'propMC');
resultGUI_MC = pln.propMC.readExternal(foldername);
