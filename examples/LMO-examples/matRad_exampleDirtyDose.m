%% Welcome to a dirty dose testing script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What is the script about?
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the script will tell you how you can use the dirty dose objectives and
% how the dose distribution will look like if you do that
% follow the steps and you will succeed :)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% open matRad files and delete everything in the workspace
matRad_rc
clear; clc;

%% load an easy phantom like the TG119.mat
load("TG119.mat")

%% Plan and Geometry
% choose your modality... but seriously choose protons!!
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

% Dirty Dose Calculation --> adding dirty dose. The number describes your
% LET threshold -> that you can change but everything else has to stay like
% this
dij = matRad_calcDirtyDose(2,dij);

%% Optimization
% yeyy only the optimization has to be done
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

% IMPORTANT: If you want to save your results, rename them! Maybe they will
% get overwritten later

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
%% Adding dirty dose
% the cst is important! 
% for adding a dirty dose objective choose a structure like the Core

% REMEMBER: Always set dirty dose objectives as a secondary objective and
% have a dose objective as your first objective
% You can change the penalty and the prescribed dirtydose if you like
cst{1,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(300,0));

%% Generate the Geometry again
stf = matRad_generateStf(ct,cst,pln);

%% Dose calculation
% Dij Calculation --> only for dose
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% Dirty Dose Calculation --> adding dirty dose. The number describes your
% LET threshold -> that you can change but everything else has to stay like
% this
dij = matRad_calcDirtyDose(2,dij);

%% Optimization
% yeyy only the optimization has to be done
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
% make sure to delete cst{1,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(300,0));
clear('cst')
load("TG119.mat")
%% Adding dirty dose
% the cst is important! 
% for adding a dirty dose objective choose a structure like the OuterTarget

% REMEMBER: Always set dirty dose objectives as a secondary objective and
% have a dose objective as your first objective
% You can change the penalty and the prescribed dirtydose if you like
cst{2,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(300,4));

%% Generate the Geometry again
stf = matRad_generateStf(ct,cst,pln);

%% Dose calculation
% Dij Calculation --> only for dose
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% Dirty Dose Calculation --> adding dirty dose. The number describes your
% LET threshold -> that you can change but everything else has to stay like
% this
dij = matRad_calcDirtyDose(2,dij);

%% Optimization
% yeyy only the optimization has to be done
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Well done your second dirty dose calculation is ready!
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
title('dirtyDose with DD objective in Target')
zoom(4)

%% Now with the Core and Target
% make sure to delete cst{2,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(300,4));
clear('cst')
load("TG119.mat")
%% Adding dirty dose
% the cst is important! 
% for adding a dirty dose objective choose two structures like the Core and OuterTarget

% REMEMBER: Always set dirty dose objectives as a secondary objective and
% have a dose objective as your first objective
% You can change the penalty and the prescribed dirtydose if you like
cst{1,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(300,0));
cst{2,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(300,4));

%% Generate the Geometry again
stf = matRad_generateStf(ct,cst,pln);

%% Dose calculation
% Dij Calculation --> only for dose
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% Dirty Dose Calculation --> adding dirty dose. The number describes your
% LET threshold -> that you can change but everything else has to stay like
% this
dij = matRad_calcDirtyDose(2,dij);

%% Optimization
% yeyy only the optimization has to be done
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Well done your third dirty dose calculation is ready!
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
title('dirtyDose with DD objective in Core and Target')
zoom(4)
