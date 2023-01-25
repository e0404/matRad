%% Example: Proton Treatment Plan with subsequent Isocenter shift
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
% (ii) how to setup a helium dose calculation 
% (iii) how to inversely optimize the pencil beam intensities directly from command window in MATLAB. 

%% set matRad runtime configuration
matRad_rc

%% Patient Data Import
load('BOXPHANTOM.mat');

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
pln.radiationMode = 'helium';        
pln.machine       = 'Generic';

%%
% for particles it is possible to also calculate the LET disutribution
% alongside the physical dose. Therefore you need to activate the
% corresponding option during dose calculcation
pln.propDoseCalc.calcLET = true;
                                       
%%
% Now we have to set the remaining plan parameters.
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [0];
pln.propStf.couchAngles   = [0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% Define the flavor of biological optimization for treatment planning along
% with the quantity that should be used for optimization. As we use helium,
% we follow a data-driven RBE parametrization to obtbain the variable 
% relative biological effectiveness. 

quantityOpt   = 'RBExD';            % either  physicalDose / effect / RBExD
modelName     = 'HEL';              % none: for photons, protons, carbon                   constRBE: constant RBE model
                                    % MCN: McNamara-variable RBE model for protons         WED: Wedenberg-variable RBE model for protons 
                                    % LEM: Local Effect Model for carbon ions              HEL: data-driven RBE parametrization for helium
% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen');

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
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot the Resulting Dose Slice
% Let's plot the transversal iso-center dose slice
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
figure
subplot(121),imagesc(resultGUI.physicalDose(:,:,slice)),colorbar,colormap(jet),title('physical dose')
subplot(122),imagesc(resultGUI.RBExD(:,:,slice)),colorbar,colormap(jet),title('RBExDose')



