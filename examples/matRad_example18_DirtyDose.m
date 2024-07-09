%% Welcome to a dirty dose example script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% open matRad files and delete everything in the workspace
matRad_rc

%% load an easy phantom like the TG119.mat
load("TG119.mat")

%% Plan and Geometry
pln.radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln.machine         = 'Generic';
pln.numOfFractions  = 30;


% beam geometry settings
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = [45 0 -45];
pln.propStf.couchAngles     = zeros(numel(pln.propStf.gantryAngles),1);  
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln.propDoseCalc.calcLET = 1;   % very important, don't forget that one!

pln.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.propOpt.spatioTemp      = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'RBExD';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,scenGenType);
stf = matRad_generateStf(ct,cst,pln);

%% Dose calculation
% Dij Calculation --> only for dose
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% Dirty Dose Calculation --> compute dirty dose influence matrix
%  arg: LET threshold (keV/um)
dij = matRad_calcDirtyDose(2,dij);

%% Optimization
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plotting
% Let's see how it looks
cube = resultGUI.physicalDose;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];
figure
subplot(2,1,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[],[]);
title('physicalDose')
zoom(4)

cube = resultGUI.dirtyDose;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];
subplot(2,1,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[],[]);
title('dirtyDose without DD objective')
zoom(4)

%% Adding dirty dose objective
% adding a dirty dose objective to the Core
cst{1,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(300,0));

%% Optimization
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Well done your first dirty dose calculation is ready!
% Let's see how it looks
cube = resultGUI.physicalDose;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];
figure
subplot(2,1,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[],[]);
title('physicalDose')
zoom(4)

cube = resultGUI.dirtyDose;
plane = 3;
slice = 80;
doseWindow = [min(cube(:)) max(cube(:))];
subplot(2,1,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[],[]);
title('dirtyDose with DD objective in Core')
zoom(4)

%% Now with the Target
%% Adding dirty dose
% for adding a dirty dose objective to OuterTarget
cst{2,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(300,20));
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Now with the Core and Target
%% Adding dirty dose
% for adding a dirty dose objective choose two structures: Core and OuterTarget
cst{1,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(300,0));
cst{2,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(300,20));
%% Optimization
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
