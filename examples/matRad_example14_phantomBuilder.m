%% Example: Generate your own phantom geometry
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team. 
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
% (i) how to create arbitrary ct data (resolution, ct numbers)
% (ii) how to create a cst structure containing the volume of interests of the phantom
% (iii) generate a treatment plan for this phantom

%% set matRad runtime configuration
clear all;
matRad_rc; %If this throws an error, run it from the parent directory first to set the paths

%% Create a CT image series
xDim = 200;
yDim = 200;
zDim = 100;

ct.cubeDim      = [yDim xDim zDim]; % second cube dimension represents the x-coordinate
ct.resolution.x = 2;
ct.resolution.y = 2;
ct.resolution.z = 3;
ct.numOfCtScen  = 1;
 
% create a ct image series with zeros - it will be filled later
ct.cubeHU{1} = ones(ct.cubeDim) * -1000;
%% Create the VOI data for the phantom
%This uses the phantombuilder class, which helps to easily implement simple
%objectives
builder = matRad_PhantomBuilder(ct);

%To add a VOI there are (so far) two different options
%either a cubic or a spherical Volume (either OAR or Target)
%NOTE: The order in which the objectives are added, matters for the order
%in the cst we will introduce later on!

%define objectives for the VOI
objective1 = struct(DoseObjectives.matRad_SquaredOverdosing(400,0));
objective2 = struct(DoseObjectives.matRad_SquaredOverdosing(10,0));
objective3 = struct(DoseObjectives.matRad_SquaredDeviation(800,45));

builder.addCubicOAR('Volume1',[60,30,60],'offset',[0 -15 0],'objectives',objective1);
builder.addCubicOAR('Volume2',[60,30,60],'offset',[0 15 0],'objectives',objective2);
builder.addSphericalTarget('Volume3',20,'objectives',objective3);

% This adds two Cubic OAR and one Spherical Target in the middle
% For Cubic Volumes a name (here Volume1 and Volume2) has to be specified,
% as well as the dimension of the Volume.
% For spherical Volumes a name has to be specified, as well as the radius
% of the sphere
%(Optional) To move the VOI from the center of the ct an offset can be set.
%(Optional) The objectives can already be set here, however this can also
%be don later on


%% Create the cst and update the ct with the generated Volumes 

builder.updatecst()

builder.updatect()

%% Get the ct and cst (stored as properties of the phantomBuilder class)
cst = builder.cst;
ct = builder.ct;
%%
display(cst);
%%
display(ct);

%if the priorities should be changed
cst{1,5}.Priority = 3;
cst{2,5}.Priority = 2;
cst{3,5}.Priority = 1,


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

%%
% Define the biological optimization model for treatment planning along
% with the quantity that should be used for optimization. Possible model values 
% are:
% 'none':     physical optimization;
% 'constRBE': constant RBE of 1.1; 
% 'MCN':      McNamara-variable RBE model for protons; 
% 'WED':      Wedenberg-variable RBE model for protons
% 'LEM':      Local Effect Model 
% and possible quantityOpt are 'physicalDose', 'effect' or 'RBExD'.
modelName    = 'none';
quantityOpt  = 'physicalDose';                                             

%%
% The remaining plan parameters are set like in the previous example files
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [0:70:355];
pln.propStf.couchAngles   = zeros(size(pln.propStf.gantryAngles));
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% retrieve nominal scenario for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen'); 

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
dij = matRad_calcPhotonDose(ct,stf,pln,cst);

%% Export dij matrix
matRad_exportDij('dij.bin',dij,stf);

%% Inverse Optimization for intensity-modulated photon therapy
% The goal of the fluence optimization is to find a set of bixel/spot 
% weights which yield the best possible dose distribution according to the
% clinical objectives and constraints underlying the radiation treatment.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot the resulting dose slice
plane      = 3;
slice      = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
doseWindow = [0 max([resultGUI.physicalDose(:)])];

figure,title('phantom plan')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],[],colorcube,[],doseWindow,[]);


%% 
% We export the the created phantom & dose as dicom. This is handled by the 
% class matRad_DicomExporter. When no arguments are given, the exporter searches
% the workspace itself for matRad-structures. The output directory can be set by
% the property dicomDir. While the different DICOM datasets (ct, RTStruct, etc) 
% can be exported individually, we call the wrapper to do all possible exports.
dcmExport = matRad_DicomExporter();
dcmExport.dicomDir = [pwd filesep 'dicomExport'];
dcmExport.matRad_exportDicom();

