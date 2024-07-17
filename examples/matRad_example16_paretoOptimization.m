%% Example: Photon Treatment Plan
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
% (ii) how to setup a photon dose calculation and 
% (iii) how to inversely optimize beamlet intensities
% (iv) how to visually and quantitatively evaluate the result
%global timed;
%timed = [];
%% Patient Data Import
% Let's begin with a clear Matlab environment. Then, import the TG119
% phantom into your workspace. The phantom is comprised of a 'ct' and 'cst'
% structure defining the CT images and the structure set. Make sure the 
% matRad root directory with all its subdirectories is added to the Matlab 
% search path.

matRad_rc; %If this throws an error, run it from the parent directory first to set the paths
matRad_rc;
matRad_cfg = MatRad_Config.instance();
%matRad_cfg.propOpt.defaultMaxIter = 50000;
%%
load('TG119.mat');
%%
cst{1,6}{1} = struct(DoseObjectives.matRad_EUD(1000,0));
cst{3,6}{1} = struct(DoseObjectives.matRad_MeanDose(1000,0));
%%
%cst{2,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(45,55));

%pln.radiationMode = 'protons';  
pln.radiationMode = 'photons';
pln.machine       = 'Generic';

quantityOpt    = 'physicalDose';                                     
modelName      = 'none';  

%pln.propDoseCalc.calcLET = 0;

pln.numOfFractions         = 30;
pln.propStf.gantryAngles   = [0:40:359];
pln.propStf.couchAngles    = zeros(1,numel(pln.propStf.gantryAngles));
pln.propStf.bixelWidth     = 5;


pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

%%
% Enable sequencing and disable direct aperture optimization (DAO) for now.
% A DAO optimization is shown in a seperate example.
pln.propOpt.runSequencing = 1;
pln.propOpt.runDAO        = 0;

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen');

%%
% and et voila our treatment plan structure is ready. Lets have a look:
display(pln);


%% Generate Beam Geometry STF
% The steering file struct comprises the complete beam geometry along with 
% ray position, pencil beam positions and energies, source to axis distance (SAD) etc.
stf = matRad_generateStf(ct,cst,pln);

%%
% Let's display the beam geometry information of the 6th beam
%display(stf(6));

%% Dose Calculation
% Let's generate dosimetric information by pre-computing dose influence 
% matrices for unit beamlet intensities. Having dose influences available 
% allows subsequent inverse optimization.
%dij = matRad_calcParticleDose(ct,stf,pln,cst);
dij = matRad_calcPhotonDose(ct,stf,pln,cst);
%% Inverse Optimization for IMRT
% The goal of the fluence optimization is to find a set of beamlet/pencil 
% beam weights which yield the best possible dose distribution according to
% the clinical objectives and constraints underlying the radiation 
% treatment. Once the optimization has finished, trigger once the GUI to 
% visualize the optimized dose cubes.
%% Paretooptimization
% The goal of this step is to define a grid of penalty values that
% are then evaluated using the matRad_paretoGeneration method
% The VOI and their respective penalties are defined in the following way
% It is possible to have more than one objective function per VOI
% penVal stores the Grid which is then passed on. penGrid contains an
% version easier to visualize, however harder to loop over
%% Add constraints
cst{1,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(0,40));
cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(0,45));
cst{2,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(45,57));

%%
[resultGUI,retStruct] = matRad_paretoOptimization(dij,cst,pln,3);

%%
matRadGUI
%%
matRadParetoGUI
% matRad_UIInterpolation(retStruct,dij,pln,ct,matRad_setOverlapPriorities(cst),retStruct.optiProb)
