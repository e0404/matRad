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

%% Example Photon Treatment Plan
% In this example we will show how to load patient data into matRad and
% how to setup a photon dose calculation and inverse optimization directly 
% from command window in MatLab.

%% Patient Data Import
% Let's begin with a clear Matlab environment. Next, import the TG119
% phantom into your workspace. The phantom is comprised of a 'ct' and 'cst' structure defining 
% the CT images and the structure set. Make sure the matRad root directy with all its
% subdirectories is added to the Matlab search path.
clc,clear,close all
load('TG119.mat');

%%
% Let's check the two variables, we have just imported. 'ct' comprises the ct cube 
% along with some meta information describing the ct cube (cube dimensions, resolution, number of CT scenarios (usefull if
% we want to plan on a 4D ct)
ct

%%
% The 'cst' cell array defines volumes of interest/structures. Each row belongs to one certain VOI, whereas each column holds
% a different information. Specifically, the second and third column show the name and the type of the structure. The fourth column depicts a linear 
% index vector depicting voxels in the CT cube that are covered by the VOI. In total, 3 structures are defined in the cst
cst

%%
% The fifth column represents meta parameters used for optimization such as the overlap priority. A lower overlap priority indicates increased importance. In contrast,
% a higher overlap priority indicatets a strcture with lower importance. The parameters alphaX and betaX depict the tissue's photon-radiosensitivity parameter
% which are required for biological treatment planning using a variable
% RBE. Let's output the meta optimization parameter of the target structure:
ixTarget = 3;
cst{ixTarget,5}

%%
% The sixth column contains optimization information such as objective and constraints which are required to calculate the objective function value. 
% Please note, that multiple objectives/constraints can be defined for individual structures. As the rectum is an OAR, we have defined and
% squared overdosing objective so that it is considered to be expensive for the optimizer delivering more than 50 Gy to the rectum. 
cst{ixTarget,6}

%% Treatment Plan
% The next step is to define your treatment plan labeld as 'pln'. This structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.

%%
% First of all, we need to define what kind of radiation modality we would
% like to use. In this case we want to use photons.
% Then, we need to define a treatment machine to correctly load the corresponding base data. Since we provide
% generic base data we set the machine to 'Genereric. By this means the the sofware will look for
% 'photons_Generic.mat' in our root directory and will use the data provided in there for dose calculation
pln.radiationMode = 'photons';
pln.machine       = 'Generic';

%%
% Define the flavour of biological optimization for treatment planning along with the quantity that should be used for
% optimizaion. As we are using photons, simply set the parameter to 'none' indicating the physical dose should be optimized.
pln.bioOptimization = 'none';

%%
% Now we set some beam parameters. We can chose multiple angles for the
% treatment and pass these to the plan as a vector and then MatRad will now
% what to do. In this case, we chose 2 parallel beams coming from opposite
% directions and we set a 5 mm bixel width.
pln.gantryAngles = [0:45:359];
pln.couchAngles  = [0 0 0 0 0 0 0 0];
pln.bixelWidth      = 5;
pln.numOfFractions  = 30;

%%
% Obtain the number of beams and voxels from the existing variables and calculate the iso-center which is per default the mass of gravity of all target voxels.
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.voxelDimensions = ct.cubeDim;
pln.isoCenter       = ones(pln.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

%%
% Enable sequencing and disable direct aperture optimization (DAO) for now.
% A DAO optimization is shown in a seperate example. The multileaf collimator leaf sequencing algorithm for intensity modulated
% beams stratifies each beam in N static segments.
pln.runSequencing = 1;
pln.runDAO        = 0;


%%
% and et voila our treatment plan is ready. Lets have a look at it:
pln

%% Generatet Beam Geometry STF
% This acronym stands for steering file and comprises the beam geomtry along with 
% the ray and pencil beam positions
stf = matRad_generateStf(ct,cst,pln);

%%
%% Let's display the beam geomtry information of the 6th beam
stf(6)

%% Dose Calculation
% Calculate dose influence matrix for unit pencil beam intensities. 
dij       = matRad_calcPhotonDose(ct,stf,pln,cst);

%% Inverse Optimizaiton for IMPT
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot the Resulting Dose Slice
% Let's plot the transversal iso-center dose slice
slice = round(pln.isoCenter(1,3)./ct.resolution.z);
figure
imagesc(resultGUI.physicalDose(:,:,slice)),colorbar

%% Now let's create another treatment plan this time leave two beams out
pln.gantryAngles = [0:45:269];
pln.couchAngles  = [0 0 0 0 0 0];
pln.numOfBeams   = numel(pln.gantryAngles);
stf              = matRad_generateStf(ct,cst,pln);
pln.isoCenter    = stf.isoCenter;
dij              = matRad_calcPhotonDose(ct,stf,pln,cst);
resultGUI_5beam  = matRad_fluenceOptimization(dij,cst,pln);

%%  Visual Comparison of results
% Let's compare the new recalculation against the optimization result.
plane = 3;
doseWindow = [0 max([resultGUI.physicalDose(:); resultGUI_5beam.physicalDose(:)])];

figure,title('original plan')
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);
figure,title('modified plan')
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_5beam.physicalDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);

%% 
% At this point we would like to see the absolute difference of the original ptimization and the 
% recalculation. 
absDiffCube = resultGUI.physicalDose-resultGUI_5beam.physicalDose;
figure,[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],[],colorcube);

%% Evalute Some Statistics


