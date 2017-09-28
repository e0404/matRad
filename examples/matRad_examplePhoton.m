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

%%
% In this example we will show 
% (i) how to load patient data into matRad
% (ii) how to setup a photon dose calculation and 
% (iii) how to inversly optimize beamlet intensities
% (iv) how to visually and quantitatively evalute the result

%% Patient Data Import
% Let's begin with a clear Matlab environment. First, import the TG119
% phantom into your workspace. The phantom is comprised of a 'ct' and 'cst' structure defining 
% the CT images and the structure set. Make sure the matRad root directory with all its
% SUBDIRECTORIES is added to the Matlab search path.
clc,clear,close all
load('TG119.mat');

%%
% Let's check the two variables, we have just imported. First, the 'ct' variable comprises the ct cube 
% along with some meta information describing properties of the ct cube (cube dimensions,
% resolution, number of CT scenarios). Please note that mutiple ct cubes
% (e.g. 4D CT) can be stored in the cell array ct.cube{}
ct

%%
% The 'cst' cell array defines volumes of interests along with information required for optimization.
% Each row belongs to one certain VOI, whereas each column defines different proprties. Specifically, the second and third column 
% show the name and the type of the structure. The tpe can be set to OAR, TARGET or IGNORED. The fourth column depicts a linear 
% index vector depicting voxels in the CT cube that are covered by the corresponding VOI. In total, 3 structures are defined in the cst
cst

%%
% The fifth column represents meta parameters used for optimization such as the overlap priority, which can be specified in double presision. 
% A lower overlap priority indicates increased importance. In contrast, a higher overlap priority indicatets a strcture with lower importance. 
% The parameters alphaX and betaX depict the tissue's photon-radiosensitivity parameter of the linear quadratic model. These parameter
% are required for biological treatment planning using a variable RBE. Let's output the meta optimization parameter of the target structure:
ixTarget = 3;
cst{ixTarget,5}

%%
% The sixth column contains optimization information such as objectives and constraints which are required to calculate the objective function value. 
% Please note, that multiple objectives/constraints can be defined for individual structures. Here, we have defined a squared deviation objective 
% making it 'expensive/costly' for the optimizer to over and underdose the target structure (both are equaly important). 
cst{ixTarget,6}

%% Treatment Plan
% The next step is to define your treatment plan labeld as 'pln'. This structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.

%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon. In this case we want to use photons.
% Then, we need to define a treatment machine to correctly load the corresponding base data. Since we provide
% generic base data we set the machine to 'Genereric. By this means matRad will look for
% 'photons_Generic.mat' in our root directory and will use the data provided in there for dose calculation
pln.radiationMode = 'photons';  
pln.machine       = 'Generic';

%%
% Define the flavour of biological optimization along with the quantity that should be used for
% optimizaion. Possible values are (none: physical optimization; const_RBExD: constant RBE of 1.1; LEMIV_effect: 
% effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose. As we are using photons, we simply set the parameter to
% 'none' thereby indicating the physical dose should be optimized.
pln.bioOptimization = 'none';    

%%
% Now we have to set some beam parameters. We can define multiple beam angles for the
% treatment and pass these to the plan as a vector. matRad will then interpret the vector as multiple beams.
% In this case, we define linear spaced beams from 0 degree to 359 degree in
% 30 degree steps. This results in 12 beams. All corresponding couch angles are set to 0 at this point.
% Moreover, we set the bixelWidth to 5, which results in a beamlet size of 5 x 5 mm in the isocenter plane. The number of fractions
% is set to 30. Be advised that matRad is always optimizing the fraction dose.
pln.gantryAngles    = [0:30:359];
pln.couchAngles     = zeros(1,numel(pln.gantryAngles));
pln.bixelWidth      = 5;
pln.numOfFractions  = 30;

%%
% Obtain the number of beams and voxels from the existing variables and calculate the iso-center which is per default
% the mass of gravity of all target voxels.
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.voxelDimensions = ct.cubeDim;
pln.isoCenter       = ones(pln.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

%%
% Enable sequencing and disable direct aperture optimization (DAO) for now.
% A DAO optimization is shown in a seperate example. The multileaf collimator leaf sequencing algorithm stratifies each beam in N static shape segments.
pln.runSequencing = 1;
pln.runDAO        = 0;

%%
% and et voila our treatment plan is ready. Lets have a look at it:
pln

%% Generatet Beam Geometry STF
% This acronym stands for steering file and comprises the complet beam geomtry along with 
% ray position, beamlet positions, source to axis distance (SAD) etc.
stf = matRad_generateStf(ct,cst,pln);

%%
% Let's display the beam geomtry information of the 6th beam
stf(6)

%% Dose Calculation
% Lets generate dosimetric information by pre-computing dose influence matrices for unit beamlet intensities. Having dose influences available allows then 
% later on inverse optimization.
dij       = matRad_calcPhotonDose(ct,stf,pln,cst);

%% Inverse Optimizaiton for IMRT
% The goal of the fluence optimization is to find a set of beamlet/pencil beam weights which yield the best possible
% dose distribution according to the clinical objectives and constraints underlying the radiation treatment.
% Once the optimization has finished, trigger once the GUI to visualize the optimized dose cubes.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
matRadGUI

%% Plot the Resulting Dose Slice
% Let's plot the transversal iso-center dose slice
slice = round(pln.isoCenter(1,3)./ct.resolution.z);
figure
imagesc(resultGUI.physicalDose(:,:,slice)),colorbar, colormap(jet)

%% Now let's create another treatment plan but this time use a coarser beam spacing.
% Instead of 30 degree spacing use a 50 degree geantry beam spacing
pln.gantryAngles = [0:50:359];
pln.couchAngles  = zeros(1,numel(pln.gantryAngles));
pln.numOfBeams   = numel(pln.gantryAngles);
stf              = matRad_generateStf(ct,cst,pln);
pln.isoCenter    = stf.isoCenter;
dij              = matRad_calcPhotonDose(ct,stf,pln,cst);
resultGUI_coarse = matRad_fluenceOptimization(dij,cst,pln);

%%  Visual Comparison of results
% Let's compare the new recalculation against the optimization result.
% Check if you have added all SUBDIRECTORIES to the MATLAB search path,
% otherwise it will not find the plotting function
plane      = 3;
doseWindow = [0 max([resultGUI.physicalDose(:); resultGUI_coarse.physicalDose(:)])];

figure,title('original plan - fine beam spacing')
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);
figure,title('modified plan - coarse beam spacing')
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_coarse.physicalDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);

%% 
% At this point we would like to see the absolute difference of the first optimization (finer beam spacing) and the 
% second optimization (coarser beam spacing)
absDiffCube = resultGUI.physicalDose-resultGUI_coarse.physicalDose;
figure,title( 'fine beam spacing plan - coarse beam spacing plan')
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],[],colorcube);

%% Obtain dose statistics
% Two more columns will be added to the cst structure depicting the DVH and
% standard dose statistics such as D95,D98, mean dose, max dose etc.
cst        = matRad_indicatorWrapper(cst,pln,resultGUI);
cst_coarse = matRad_indicatorWrapper(cst,pln,resultGUI_coarse);

%%
% The treatment plan using more beams should in principle result in a
% better OAR sparing. Therefore lets have a look at the D95 of the OAR of both plans
ixOAR = 2;
cst{ixOAR,9}{1}.D95
cst_coarse{ixOAR,9}{1}.D95

%%
% Indeed, the coarse beam plan yields a higher D95 in the OAR. Finally, lets plot the DVH 
matRad_showDVH(cst,pln)



