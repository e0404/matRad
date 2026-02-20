% matRad script
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set matRad runtime configuration
matRad_rc

%% load patient data, i.e. ct, voi, cst
load TG119.mat
%load HEAD_AND_NECK
%load PROSTATE.mat
%load LIVER.mat
%load BOXPHANTOM.mat

% meta information for treatment plan
pln.numOfFractions  = 30;
pln.radiationMode   = 'photons';            % either photons / protons / helium / carbon / brachy / VHEE
pln.machine         = 'Generic';            % generic for RT / LDR or HDR for BT / Generic or Focused for VHEE

pln.bioModel = 'none';      % none: for all                                 % constRBE: constant RBE for photons and protons 
                            % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                            % LEM: Local Effect Model for carbon ions       % HEL: data-driven RBE parametrization for helium

pln.multScen = 'nomScen';   % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'         

% beam geometry settings
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = [0:72:359]; % [°] ;
pln.propStf.couchAngles     = [0 0 0 0 0]; % [°] ; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

% optimization settings
pln.propOpt.quantityOpt     = 'physicalDose';   % Quantity to optimizer (could also be RBExDose, BED, effect)
pln.propOpt.optimizer       = 'IPOPT';          % We can also utilize 'fmincon' from Matlab's optimization toolbox
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propSeq.runSequencing   = true;  % true: run sequencing, false: don't / will be ignored for particles and also triggered by runDAO below


%% initial visualization and change objective function settings if desired
matRadGUI

%% generate steering file 
stf = matRad_generateStf(ct,cst,pln);

%% dose calculation
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

%% inverse planning for imrt
resultGUI  = matRad_fluenceOptimization(dij,cst,pln);

%% sequencing
resultGUI = matRad_sequencing(resultGUI,stf,dij,pln);


%% DAO
if strcmp(pln.radiationMode,'photons') && pln.propOpt.runDAO
   resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);
   matRad_visApertureInfo(resultGUI.apertureInfo);
end

%% start gui for visualization of result
matRadGUI

%% indicator calculation and show DVH and QI
resultGUI = matRad_planAnalysis(resultGUI,ct,cst,stf,pln);

