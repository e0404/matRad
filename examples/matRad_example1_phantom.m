%% Example: Generate your own phantom geometry
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In this example we will show
% (i) how to create arbitrary ct data (resolution, ct numbers)
% (ii) how to create a cst structure containing the volume of interests of the phantom
% (iii) generate a treatment plan for this phantom

%% set matRad runtime configuration
%clear all; %somewhat needed for the phantom builder
matRad_rc; %If this throws an error, run it from the parent directory first to set the paths
matRad_cfg = MatRad_Config.instance(); %This creates a matRad-Config object holding global configuration parameters

%% Create a CT image series

ctDim = [200,200,100]; % x,y,z dimensions
ctResolution = [2,2,3]; % x,y,z the same here!

%This uses the phantombuilder class, which helps to easily implement a 
%water phantom containing geometrical 3D shapes as targets and organs
builder = matRad_PhantomBuilder(ctDim,ctResolution,1);

%% Create the VOI data for the phantom
%To add a VOI there are (so far) two different options
%either a Box or a spherical Volume (either OAR or Target)
%NOTE: The order in which the objectives are initialized matters!
% In case of overlaps in the objectives, the firstly created objectives have
% a higher priority! This means that if two VOI have an overlap with
% different HU, then the value of the firstly initialized objective will be
% set in the overlap region


%define objectives for the VOI

objective1 = struct(DoseObjectives.matRad_SquaredDeviation(800,45));
objective2 = struct(DoseObjectives.matRad_SquaredOverdosing(400,0));
objective3 = struct(DoseObjectives.matRad_SquaredOverdosing(10,0));     

builder.addSphericalTarget('Volume1',20,'objectives',objective1,'HU',0);
builder.addBoxOAR('Volume2',[60,30,60],'offset',[0 -15 0],'objectives',objective2);
builder.addBoxOAR('Volume3',[60,30,60],'offset',[0 15 0],'objectives',objective3);

% This adds two Box OAR and one Spherical Target in the middle
% For Box Volumes a name (here Volume2 and Volume3) has to be specified,
% as well as the dimension of the Volume.
% For spherical Volumes a name has to be specified, as well as the radius
% of the sphere
%(Optional) To move the VOI from the center of the ct an offset can be set.
%(Optional) The objectives can already be set here, however this can also
%be done later on
%(Optional) The HU of the VOI can be set (normal value: 0 (water)) 


%% Get the ct and cst (stored as properties of the phantomBuilder class)

[ct,cst] = builder.getctcst();

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines the most
% important cornerstones of your treatment plan.
%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon. In this
% example we would like to use photons for treatment planning. Next, we
% need to define a treatment machine to correctly load the corresponding 
% base data. matRad features generic base data in the file
% 'photons_Generic.mat'; consequently the machine has to be set to 'Generic'
pln.radiationMode = 'photons';            
pln.machine       = 'Generic';

% Define the biological optimization model for treatment planning.
% Examples for possible models are:
% 'none':     no biological model (physical dose only);
% 'constRBE': constant RBE of 1.1; 
% 'MCN':      McNamara-variable RBE model for protons; 
% 'WED':      Wedenberg-variable RBE model for protons
% 'LEM':      Local Effect Model (based on precomputed kernels)
pln.bioModel = 'none';

% matRad allows for an uncertainty scenario model to be respected in planning, 
% which compute a number of error scenarios which can be used in optimization.
% Examples for possible scenario models
% 'nomScen':    Only Nominal Scenario is considered
% 'wcScen':     Worst Case Scenarios are generated
% 'impScen':    Importance weighted gridded scenarios are generated
% 'rndScen':    Randomly sampled scenarios are generated
pln.multScen = 'nomScen';

% Number of fractions
% Matlab considers prescriptions for the total plan, but shows fraction dose
pln.numOfFractions        = 30;

% Settings for the irradiation geometry
pln.propStf.gantryAngles  = [0:70:355];
pln.propStf.couchAngles   = zeros(size(pln.propStf.gantryAngles));
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

% Settings for Optimization
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

%% Visualization in GUI
matRadGUI;

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% Inverse Optimization for intensity-modulated photon therapy
% The goal of the fluence optimization is to find a set of bixel/spot 
% weights which yield the best possible dose distribution according to the
% clinical objectives and constraints underlying the radiation treatment.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Visualization
matRadGUI
%% Plot the resulting dose slice
plane      = 3;
slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
slice = slice(3);
doseWindow = [0 max([resultGUI.physicalDose(:)])];

figure,title('phantom plan')
matRad_plotSlice(ct, 'axesHandle', gca, 'cst', cst, 'cubeIdx', 1, 'dose', resultGUI.physicalDose, 'plane', plane, 'slice', slice, 'contourColorMap', colorcube, 'doseWindow', doseWindow);
%matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],[],colorcube,[],doseWindow,[]);

%% 
% We export the the created phantom & dose as dicom. This is handled by the 
% class matRad_DicomExporter. When no arguments are given, the exporter searches
% the workspace itself for matRad-structures. The output directory can be set by
% the property dicomDir. While the different DICOM datasets (ct, RTStruct, etc) 
% can be exported individually, we call the wrapper to do all possible exports.
dcmExport = matRad_DicomExporter();
dcmExport.dicomDir = [matRad_cfg.primaryUserFolder filesep 'dicomExport'];
dcmExport.matRad_exportDicom();

