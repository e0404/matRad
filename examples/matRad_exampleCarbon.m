%% Example: Carbon Ion Treatment Plan
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

%%
% In this example we will demonstrate:
% (i)  how to load one of the open source patient datasets into matRad 
% (ii) how to setup a carbon ion treatment plan including variable RBE optimization
% (iii) perform a recalculation assuming a different tissue radio-sensitivity

%% Import
% Let's begin with a clear Matlab environment. Next, import the liver patient data
% into your workspace. The patient is comprised of a 'ct' and 'cst' structure defining 
% the CT images and the structure set. Make sure the matRad root directy with all its
% subdirectories is added to the Matlab search path.
clc,clear,close all
load('LIVER.mat');

%%
% Let's check the two variables, we have just imported. First, the 'ct' variable comprises the ct cube 
% along with some meta information describing properties of the ct cube (cube dimensions,
% resolution, number of CT scenarios). Please note that mutiple ct cubes
% (e.g. 4D CT) can be stored in the cell array ct.cube{}
ct

%%
% The 'cst' cell array defines volumes of interest along with information required for optimization.
% Each row belongs to one certain VOI, whereas each column defines different proprties. Specifically, the second and third column 
% show the name and the type of the structure. The tpe can be set to OAR, TARGET or IGNORED. The fourth column depicts a linear 
% index vector depicting voxels in the CT cube that are covered by the VOI. In total, 17 structures are defined in the cst
cst
%%
% The fifth column represents meta parameters used for optimization such as the overlap priority, which can be specified in double presision. A lower overlap priority indicates increased importance. In contrast,
% a higher overlap priority indicatets a strcture with lower importance. The parameters alphaX and betaX depict the tissue's photon-radiosensitivity parameter
% which are required for biological treatment planning using a variable RBE. Since the PTV, is stored in the tenth column, lets output the PTVs meta information
cst{10,5}

%%
% The sixth column contains optimization information such as objectives and constraints which are required to calculate the objective function value. 
% Please note, that multiple objectives/constraints can be defined for individual structures. Here, we have defined a squared deviation
% objective for the target thereby equaly penalizing over and underdosage.
cst{10,6}

%%
% Now lets increase the importance of the square deviation objective to 2000
cst{10,6}.penalty = 2000;

%% Treatment Plan
% The next step is to define your treatment plan labeld as 'pln'. This structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.
%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon. In this
% example we would like to use carbon ions for treatment planning.
% Then, we need to define a treatment machine to correctly load the corresponding base data. Since we provide
% generic base data we set the machine to 'Genereric. By this means matRad will look for
% 'carbon_Generic.mat' in our root directory and will use the data provided in there for dose calculation
pln.radiationMode = 'carbon';            
pln.machine       = 'Generic';

%%
% Define the flavour of biological optimization for treatment planning along with the quantity that should be used for
% optimizaion. Possible values are (none: physical optimization; const_RBExD: constant RBE of 1.1; LEMIV_effect: 
% effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose. As we use carbon ions, we what to use 
% base data from the local effect model IV and want to optimize the RBE-weighted dose. Therefore we set bioOptimization to LEMIV_RBExD
pln.bioOptimization = 'LEMIV_RBExD';                                              

%%
% Now we have to set some beam parameters. We can define multiple beam angles for the
% treatment and pass these to the plan as a vector. matRad will then interpret the vector as multiple beams. However, now we define a single beam
% with a gantry angle of 315 degree and set the corresponding couch angle to 0 degree. Furthermore, we want the lateral pencil beam spacing in x and y to be 3 mm in the iso-center.
% In total we are using 30 fractions. It is noteworthy that matRad is always optimizing the fraction dose.
pln.gantryAngles    = 315;
pln.couchAngles     = 0;
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
% and schwuppdiwupp - our treatment plan is ready. Lets have a look at it:
pln

%% Generatet Beam Geometry STF
% This acronym stands for steering file and comprises the beam geomtry along with 
% the ray and pencil beam positions
stf = matRad_generateStf(ct,cst,pln);

%%
% Let's display the beam geomtry information
stf

%% Dose Calculation
% Calculate dose influence matrix for unit pencil beam intensities. 
dij = matRad_calcParticleDose(ct,stf,pln,cst);

%% Inverse Optimizaiton for IMPT
% The goal of the fluence optimization is to find a set of bixel/spot weights which yield the best possible
% dose distribution according to the clinical objectives and constraints underlying the radiation treatment.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot the Resulting Dose Slice
% Let's plot the transversal iso-center dose slice
slice = round(pln.isoCenter(3)./ct.resolution.z);
figure,
imagesc(resultGUI.RBExDose (:,:,slice)),colorbar, colormap(jet)

%% Manipulate the cst
% The previous treatment plan was optimized using an photon alpha-beta ratio of 2. Now, Let's
% change the radiosensitivity by adapting alphaX. This will change the photon alpha-beta ratio
% from 2 to 10.
for i = 1:size(cst,1)
    cst{i,5}.alphaX      = 0.5;
    cst{i,5}.TissueClass = 2;
end

%% Recalculate Plan
% Let's use the existing optimized pencil beam weights and recalculate the RBE weighted dose
resultGUI_tissue = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);


%% Result Comparison
% Let's compare the new recalculation against the optimization result.
plane = 3;
doseWindow = [0 max([resultGUI.RBExDose(:); resultGUI_tissue.RBExDose(:)])];

figure,title('original plan')
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.RBExDose,plane,slice,[],[],colorcube,[],doseWindow,[]);
figure,title('manipulated plan')
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_tissue.RBExDose,plane,slice,[],[],colorcube,[],doseWindow,[]);

%% 
% At this point we would like to see the absolute difference of the original ptimization and the 
% recalculation. 
absDiffCube = resultGUI.RBExDose-resultGUI_tissue.RBExDose;
figure,title('absolute difference')
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],[],colorcube);




