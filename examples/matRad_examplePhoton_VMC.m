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
% (ii) how to setup a photon dose calculation based on the VMC++ Monte Carlo algorithm 
% (iii) how to inversely optimize the beamlet intensities directly from command window in MATLAB. 
% (iv) how to visualize the result

%% Patient Data Import
% Let's begin with a clear Matlab environment. First, import the boxphantom
% into your workspace. The phantom is comprised of a 'ct' and 'cst' structure defining 
% the CT images and the structure set. Make sure the matRad root directory with all its
% SUBDIRECTORIES is added to the Matlab search path.
clc,clear,close all
load('BOXPHANTOM.mat');

%%
% Let's check the two variables, we have just imported. First, the 'ct' variable comprises the ct cube 
% along with some meta information describing properties of the ct cube (cube dimensions,
% resolution, number of CT scenarios). Please note that multiple ct cubes
% (e.g. 4D CT) can be stored in the cell array ct.cube{}
ct

%%
% The 'cst' cell array defines volumes of interests along with information required for optimization.
% Each row belongs to one certain VOI, whereas each column defines different properties. Specifically, the second and third column 
% show the name and the type of the structure. The type can be set to OAR, TARGET or IGNORED. The fourth column depicts a linear 
% index vector depicting voxels in the CT cube that are covered by the corresponding VOI. In total, 2 structures are defined in the cst
cst

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This structure requires input from the treatment planner and defines 
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
% Define the flavor of biological optimization along with the quantity that should be used for
% optimization. Possible values are (none: physical optimization; const_RBExD: constant RBE of 1.1; LEMIV_effect: 
% effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose. As we are using photons, we simply set the parameter to
% 'none' thereby indicating the physical dose should be optimized.
pln.bioOptimization = 'none';    


%%
% Now we have to set some beam parameters. We can define multiple beam angles for the
% treatment and pass these to the plan as a vector. matRad will then interpret the vector as multiple beams.
% In this case, we define one single beam from 0 degree gantry angle and 0 degree couch angle.
% Moreover, we set the bixelWidth to 10, which results in a beamlet size of 10 x 10 mm. The number of fractions
% is set to 30. Be advised that matRad is always optimizing the fraction dose.
pln.gantryAngles    = [0];
pln.couchAngles     = [0];
pln.bixelWidth      = 10;
pln.numOfFractions  = 30;

%%
% Obtain the number of beams and voxels from the existing variables and calculate the iso-center which is per default
% the mass of gravity of all target voxels.
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.voxelDimensions = ct.cubeDim;
pln.isoCenter       = ones(pln.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

%%
% Disable sequencing and direct aperture optimization (DAO) for now.
% The application of the sequencing algorithm and DAO optimization is shown in a separate example.
pln.runSequencing = 0;
pln.runDAO        = 0;

%%
% and et voila our treatment plan is ready. Lets have a look at it:
pln

%% Generate Beam Geometry STF
% This acronym stands for steering file and comprises the complete beam geometry along with 
% ray position, beamlet positions, source to axis distance (SAD) etc.
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
% Calculate dose influence matrix for unit pencil beam intensities. using the
% VMC++ monte carlo algorithm. We define the number of photons simulated
% per beamlet to be 700. Make sure to download VMC++ files from http://www.cerr.info/download.php
% into the matRadrootDirectory\vmc++ to run the following function.
nCasePerBixel = 700;              
dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst,700,1);

%% Inverse Optimization for IMRT
% The goal of the fluence optimization is to find a set of bixel/spot weights which yield the best possible
% dose distribution according to the clinical objectives and constraints underlying the radiation treatment.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot the Resulting Dose Slice
% Just let's plot the transversal iso-center dose slice
slice = round(pln.isoCenter(1,3)./ct.resolution.z);
figure,
imagesc(resultGUI.physicalDose(:,:,slice)),colorbar, colormap(jet)

%%
% Exemplary, we show how to obtain the dose in the target and plot the histogram
ixTarget     = cst{2,4}{1};
doseInTarget = resultGUI.physicalDose(ixTarget);
figure,histogram(doseInTarget),title('dose in target'),xlabel('[Gy]')

%% Start the GUI for Visualization
matRadGUI


