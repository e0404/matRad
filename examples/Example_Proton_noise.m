% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad proton code script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Example Proton Treatment Plan with Manipulated CT values
% In this example we will show how to load patient data into matRad and
% how to setup a proton dose calculation. Next we inverslye optimize the treatment plan 
% directly from command window in MatLab. Next, we add noise to the CT cube
% and will perform a forward dose calculation. 


%% Patient Data Import
% Let's begin with a clear Matlab environment. Next, import the protstate
% case into your workspace. The patient is comprised of a 'ct' and 'cst' structure defining 
% the CT images and the structure set. Make sure the matRad root directy with all its
% subdirectories is added to the Matlab search path.
clc,clear,close all
load('PROSTATE.mat');

%%
% Let's check the two variables, we have just imported. 'ct' comprises the ct cube 
% along with some meta information describing the ct cube (cube dimensions, resolution, number of CT scenarios (usefull if
% we want to plan on a 4D ct)
ct

%%
% The 'cst' cell array defines volumes of interest/structures. Each row belongs to one certain VOI, whereas each column holds
% a different information. Specifically, the second and third column show the name and the type of the structure. The fourth column depicts a linear 
% index vector depicting voxels in the CT cube that are covered by the VOI. In total, 17 structures are defined in the cst
cst

%%
% The fifth column represents meta parameters used for optimization such as the overlap priority. A lower overlap priority indicates increased importance. In contrast,
% a higher overlap priority indicatets a strcture with lower importance. The parameters alphaX and betaX depict the tissue's photon-radiosensitivity parameter
% which are required for biological treatment planning using a variable
% RBE. Let's output the meta optimization parameter of the rectum the rectum:
cst{8,2} 
cst{8,5}

%%
% The sixth column contains optimization information such as objective and constraints which are required to calculate the objective function value. 
% Please note, that multiple objectives/constraints can be defined for individual structures. As the rectum is an OAR, we have defined and
% squared overdosing objective so that it is considered to be expensive for the optimizer delivering more than 50 Gy to the rectum. 
cst{8,6}


%% Treatment Plan
% The next step is to define your treatment plan labeld as 'pln'. This structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.

%%
% First of all, we need to define what kind of radiation modality we would like to use.
% Then, we need to define a treatment machine to correctly load the corresponding base data. Since we provide
% generic base data we set the machine to 'Genereric. By this means the the sofware will look for
% 'proton_Generic.mat' in our root directory and will use the data provided in there for dose calculation
pln.radiationMode = 'protons';
pln.machine       = 'Generic';

%%
% Define the flavour of biological optimization for treatment planning along with the quantity that should be used for
% optimizaion. const_RBExD results in using a constant RBE of 1.1 and use the RBE-weighted dose during optimization.
pln.bioOptimization = 'const_RBExD';

%%
% Now we have to further specify the treatment plan. We define two opposing beams. For the first beam we set the gantry angle to 90 degree and the 
% corresponding couch angle to 0 degree. The second beam possess a gantry angle of 270 degree and a couch angle of 0 degree. 
% Furthermore, we want the lateral pencil beam spacing in x and y to be 3 mm in the iso-center slice.
% In total we are using 30 fractions. It is noteworthy that matRad is always optimizing the fraction dose.
pln.gantryAngles    = [90 270];
pln.couchAngles     = [0 0];
pln.bixelWidth      = 3;
pln.numOfFractions  = 30;

%%
% Obtain the number of beams and voxels from the existing variables and calculate the iso-center which is per default the mass of gravity of all target voxels.
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.voxelDimensions = ct.cubeDim;
pln.isoCenter       = ones(pln.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

%%
% Disable sequencing and direct aperture optimization, as we have a particle plan.
pln.runDAO        = 0;
pln.runSequencing = 0;

%%
% and et voila our treatment plan is ready. Lets have a look at it:
pln

%% Generatet Beam Geometry STF
% This acronym stands for steering file and comprises the beam geomtry along with 
% the ray and pencil beam positions
stf = matRad_generateStf(ct,cst,pln);

%% Let's display the beam geomtry information of the second beam
stf(2)

%% Dose Calculation
% Calculate dose influence matrix for unit pencil beam intensities. 
dij = matRad_calcParticleDose(ct,stf,pln,cst);

%% Inverse Optimizaiton for IMPT
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot the Resulting Dose Slice
% Let's plot the transversal iso-center dose slice
slice = round(pln.isoCenter(1,3)./ct.resolution.z);
figure
imagesc(resultGUI.RBExDose(:,:,slice)),colorbar

%%
% Now let's simulate an range overshoot by scaling the relative stopping power cube 3.5% percent
ct_manip         = ct;
noise            = ct.cube{1}  .* 0.035;
ct_manip.cube{1} = ct_manip.cube{1} - noise;

%% Recalculate Plan
% Let's use the existing optimized pencil beam weights and recalculate the RBE weighted dose
resultGUI_noise = matRad_calcDoseDirect(ct_manip,stf,pln,cst,resultGUI.w);

%%  Visual Comparison of results
% Let's compare the new recalculation against the optimization result.
plane = 3;
doseWindow = [0 max([resultGUI.RBExDose(:); resultGUI_noise.RBExDose(:)])];

figure,title('original plan')
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.RBExDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);
figure,title('manipulated plan')
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct_manip,cst,1,resultGUI_noise.RBExDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);


%% Quantitative Comparison of results
% Compare the two dose cubes using a gamma-index analysis. To do so, we need
% to define thresholds for passing the gamma-index test. In this example,
% we set dose difference to 2 % and the distance-to criteria to 2 mm.
% Additionaly, we set the number of interpolations for the gamma index
% calculation to 1. 
criteria = [2, 2];
n = 1;

[gammaCube,gammaPassRateCell] = matRad_gammaIndex(...
    resultGUI_noise.RBExDose,resultGUI.RBExDose,...
    [ct.resolution.x, ct.resolution.y, ct.resolution.z],...
    slice,criteria,n,'global',cst);




