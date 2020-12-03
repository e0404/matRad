%% Example: Carbon Ion Treatment Plan
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
% (ii) how to setup a carbon ion dose calculation plan including variable RBE optimization
% (iii) how to inversely optimize the pencil beam intensities based on the
% RBE-weighted dose
% (iv) how to inversely optimize the pencil beam intensities based on the
% biological effect
% (v) how to change the tissues' radiobiological characteristics
% (vi) how to recalculated the dose considering the previously optimized pencil beam intensities
% (vii) how to compare the two results

%% Patient Data Import
% Let's begin with a clear Matlab environment and import the liver
% patient into your workspace.

matRad_rc; %If this throws an error, run it from the parent directory first to set the paths

load('LIVER.mat');

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines the most
% important cornerstones of your treatment plan.
%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon. In this
% example we would like to use carbon ions for treatment planning. Next, we
% need to define a treatment machine to correctly load the corresponding 
% base data. matRad features generic base data in the file
% 'carbon_Generic.mat'; consequently the machine has to be set accordingly
pln.radiationMode = 'carbon';            
pln.machine       = 'Generic';

%%
% Define the flavor of biological optimization for treatment planning along
% with the quantity that should be used for optimization. Possible values 
% are (none: physical optimization; const_RBExD: constant RBE of 1.1; 
% LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of 
% RBE-weighted dose. As we use carbon ions, we decide to use base data from 
% the local effect model IV and want to optimize the RBE-weighted dose. 
% Therefore we set bioOptimization to LEMIV_RBExD
pln.propOpt.bioOptimization = 'LEMIV_RBExD';                                              

%%
% The remaining plan parameters are set like in the previous example files
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = 315;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 3;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%%
% Let's have a closer look on the stf.ray sub-structure which contains the 
% actual  beam/ray geometry information. For illustration purposes we want 
% to show the last ray. Besides geometrical  information about the position 
% and orientation of the ray, we can also find pencil beam information. If 
% the ray coincides with the target, pencil beams were defined along the 
% ray from target entry to target exit. 
display(stf.ray(end));

%%
% Here are the energies selected on the last ray: 
display(stf.ray(end).energy);

%% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

%% Inverse Optimization  for IMPT based on RBE-weighted dose
% The goal of the fluence optimization is to find a set of bixel/spot 
% weights which yield the best possible dose distribution according to the
% clinical objectives and constraints underlying the radiation treatment.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot the Resulting Dose Slice
% Let's plot the transversal iso-center dose slice
slice = round(pln.propStf.isoCenter(3)./ct.resolution.z);
figure,
imagesc(resultGUI.RBExDose (:,:,slice)),colorbar, colormap(jet);

%% Inverse Optimization  for IMPT based on biological effect
% To perform a dose optimization for carbon ions we can also use the
% biological effect instead of the RBE-weighted dose. Therefore we have to
% change the optimization mode and restart the optimization
pln.propOpt.bioOptimization = 'LEMIV_effect'; 
resultGUI_effect = matRad_fluenceOptimization(dij,cst,pln);

%% Visualize differences
% Through optimzation based on the biological effect we obtain a slightly
% different dose distribution as visualized by the following dose
% difference map
figure;
imagesc(resultGUI.RBExDose (:,:,slice)-resultGUI_effect.RBExDose(:,:,slice));
colorbar;
colormap(jet);

%% Change Radiosensitivity
% The previous treatment plan was optimized using an photon alpha-beta 
% ratio of 2 for all tissues. Now, Let's change the radiosensitivity by 
% adapting alphaX. This will change the photon alpha-beta ratio
% from 2 to 10.
for i = 1:size(cst,1)
    cst{i,5}.alphaX      = 0.5;
    cst{i,5}.TissueClass = 2;
end

%% Recalculate Plan
% Let's use the existing optimized pencil beam weights and recalculate the RBE weighted dose
resultGUI_tissue = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);

%% Result Comparison
% Let's compare the new recalculation against the optimization result.
plane = 3;
doseWindow = [0 max([resultGUI_effect.RBExDose(:); resultGUI_tissue.RBExDose(:)])];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_effect.RBExDose,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('original plan')
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_tissue.RBExDose,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
%% 
% At this point we would like to see the absolute difference of the original optimization and the 
% recalculation. 
absDiffCube = resultGUI_effect.RBExDose-resultGUI_tissue.RBExDose;
figure,
matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],[],colorcube);
title('absolute difference')
%%
% Plot both doses with absolute difference and gamma analysis
[gammaCube,gammaPassRate,hfigure]=matRad_compareDose(resultGUI_effect.RBExDose, resultGUI.RBExDose, ct, cst,[1 1 1],'on');



