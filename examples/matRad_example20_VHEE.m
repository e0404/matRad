%% Example: VHEE Treatment Plan
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this example contributors 
% Authors : F. D'Andrea ; A. Bennan ; L. Ermeneux ; N. Wahl 
%
% Based on implementation proposed by M. Sitarz et al. (doi:10.1002/mp.17392)
% Generic machine using FermiEyges model based on work of M.G. Ronga et al. (doi:10.1002/mp.16697)
% Applied matRad for a VHEE study as described by F. D'andrea et al. (doi:10.1016/j.phro.2025.100732)
% Focused machine based on work of L. Whitmore et al. (doi:10.1038/s41598-021-93276-8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In this example we will show 
% (i) how to load patient data into matRad
% (ii) how to setup a VHEE dose calculation 
% (iii) how to inversely optimize the pencil beam intensities directly from 
%       command window in MATLAB. 

%% set matRad runtime configuration
matRad_rc; %If this throws an error, run it from the parent directory first to set the paths

%% Patient Data Import
% Let's begin with a clear Matlab environment and import the prostate
% patient into your workspace
load('PROSTATE.mat');

%% Treatment Plan
% Here, we would like to use VHEE for treatment planning. Next, we need to
% define a treatment machine to correctly load the corresponding base data.
% matRad features two base data for VHEE, a divergent beam 
% (VHEE_Generic.mat) based on a FermiEyges model that has to be called 
% through 'Generic' and a Focused beam (VHEE_Focused.mat), to be called by 
% 'Focused'.
pln.radiationMode   = 'VHEE';    % either photons / protons / helium / carbon / brachy / VHEE
pln.machine         = 'Generic'; %  Generic / Focused VHEE - (Focused still in development) 
pln.bioModel        = 'none'; % 'none' for VHEE
                                       
%% plan parameters
% Now we have to set the remaining plan parameters.
% beam geometry settings
pln.numOfFractions          = 30;
pln.propStf.energy          = 200; % set VHEE beam energy in MeV [100,150 or 200 MeV]
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = [35, 110, 180, 250, 325]; % [°] ;
pln.propStf.couchAngles     = [0 0 0 0 0]; % [°] ; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]
pln.propDoseCalc.engine = 'HongPB';

% optimization settings
pln.propOpt.quantityOpt     = 'physicalDose';   % Quantity to optimizer (could also be RBExDose, BED, effect)
pln.propOpt.optimizer       = 'IPOPT';          % We can also utilize 'fmincon' from Matlab's optimization toolbox
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propSeq.runSequencing   = false;  % true: run sequencing, false: don't / will be ignored for particles and also triggered by runDAO below

%% generate steering file 
stf = matRad_generateStf(ct,cst,pln);

%% dose calculation
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

%% inverse planning for imrt
resultGUI  = matRad_fluenceOptimization(dij,cst,pln);  % Future work - remove low weighted spots to aid MC

%% use the GUI widgets directly to visualize the result
viewer = matRad_ViewingWidget();
viewer.doseOpacity = 0.35; %lets change the doseOpacity

dvhwidget = matRad_DVHStatsWidget();
dvhwidget.selectedDisplayOption = 'physicalDose';

%% Export parameter files for a TOPAS recalculation
% set number of histories lower than default for this example (default: 1e8)
pln.propDoseCalc.numHistoriesDirect = 5e6;
pln.propDoseCalc.engine = 'TOPAS';
pln.propDoseCalc.externalCalculation = 'write';
resultGUI_MC = matRad_calcDoseForward(ct,cst,stf,pln,resultGUI.w);

