%% Example Proton Treatment Plan with Manipulated CT values
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
% (iii) how to inversly optimize directly from command window in MatLab.
% (iv) how to apply a sequencing algorithm
% (v) how to run a direct aperture optimization
% (iv) how to visually and quantitatively evalute the result

%% Patient Data Import
% Let's begin with a clear Matlab environment. First, import the head &
% neck patient into your workspace. The phantom is comprised of a 'ct' and 'cst' structure defining 
% the CT images and the structure set. Make sure the matRad root directory with all its
% SUBDIRECTORIES is added to the Matlab search path.
clc,clear,close all
load('HEAD_AND_NECK.mat');

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
% index vector depicting voxels in the CT cube that are covered by the corresponding VOI. In total, 24 structures are defined in the cst
cst

%% Treatment Plan
% The next step is to define your treatment plan labeld as 'pln'. This structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.

%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon. In this case we want to use photons.
% Then, we need to define a treatment machine to correctly load the corresponding base data. Since we provide
% generic base data we set the machine to 'Genereric. By this means matRad will look for
% 'proton_Generic.mat' in our root directory and will use the data provided in there for dose calculation
pln.radiationMode = 'photons';   % either photons / protons / carbon
pln.machine       = 'Generic';

%%
% Define the flavour of biological optimization for treatment planning along with the quantity that should be used for
% optimizaion. Possible values are (none: physical optimization; const_RBExD: constant RBE of 1.1; LEMIV_effect: 
% effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose. As we are using photons, simply set the parameter to
% 'none' thereby indicating the physical dose should be optimized.
pln.bioOptimization = 'none';    

%%
% Now we have to set some beam parameters. We can define multiple beam angles for the
% treatment and pass these to the plan as a vector.matRad will then interpret the vector as multiple beams.
% In this case, we define linear spaced beams from 0 degree to 359 degree in
% 30 degree steps. This results in 12 beams. All corresponding couch angles are set to 0 at this point.
% Moreover, we set the bixelWidth to 5, which results in a beamlet size of 5 x 5 mm. The number of fractions
% is set to 30. Be advised that matRad is always optimizing the fraction dose.
pln.gantryAngles    = [0:40:359];
pln.couchAngles     = [0 0 0 0 0 0 0 0 0];
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
% Enable sequencing and direct aperture optimization (DAO).
pln.runSequencing = 1;
pln.runDAO        = 1;

%%
% y listo our treatment plan is ready. Lets have a look at it:
pln

%% Generatet Beam Geometry STF
% This acronym stands for steering file and comprises the complet beam geomtry along with 
% ray position, beamlet positions, source to axis distance (SAD) etc.
stf = matRad_generateStf(ct,cst,pln);

%% Start the GUI for Visualization
% Show the ct cube along with contours
matRadGUI

%% Dose Calculation
% Lets generate dosimetric information by pre-computing dose influence matrices for unit beamlet intensities. Having dose influences available allows then 
% later on inverse optimization.
dij = matRad_calcPhotonDose(ct,stf,pln,cst);

%% Inverse Planning for IMRT
% The goal of the fluence optimization is to find a set of beamlet weights which yield the best possible
% dose distribution according to the predefined clinical objectives and constraints underlying the radiation treatment.
% Once the optimization has finished, trigger once the GUI to visualize the optimized dose cubes.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
matRadGUI

%% Sequencing
% This is a multileaf collimator leaf sequencing algorithm that is used in 
% order to modulate the intensity of the beams with multiple static 
% segments, so that translates each intensity map into a set of deliverable 
% aperture shapes; according to Siochi (1999).
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,5);

%% DAO - Direct Aperture Optimization
% The Direct Aperture Optimization is an automated planning system, only
% possibble for photons in which we bypass the traditional intensity
% optimization, and instead directly optimize the shapes and the weights of
% the apertures. This technique allows the user to specify the maximum
% number of apertures per beam direction, and hence provides significant
% control over the complexity of the treatment delivery.
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);
matRad_visApertureInfo(resultGUI.apertureInfo);

%% Indicator Calculation
% Here we call the indicator wrapper function that is just a wrapper
% function to calculate multiple quality indicators like D98, D95. One
% basically calculates default quality indicators for a given dose distribution.
% If you don't see clearly the images, go in the folder 'html' in your
% working folder and you should look at the images, starting with the same
% name of this file, numbered from 3 to 11.
cst = matRad_indicatorWrapper(cst,pln,resultGUI);

%% Show DVH and QI
% This last function allow you to plot the Dose Volume Histogram and the
% dose statistics of all the VOIs in side the ct.
matRad_showDVH(cst,pln)

