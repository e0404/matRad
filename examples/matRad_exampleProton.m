%% Example: Proton Treatment Plan with subsequent Isocenter shift
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
% (ii) how to setup a proton dose calculation 
% (iii) how to inversly optimize the pencil beam intensities directly from command window in MatLab. 
% (iv) how to simulate a lateral patient displacement by shifting the iso-center 
% (v) how to recalculated the dose considering the shifted geometry and the previously optimzed pencil beam intensities
% (vi) how to compare the two results

%% Patient Data Import
% Let's begin with a clear Matlab environment. First, import the prostate
% patient into your workspace. The phantom is comprised of a 'ct' and 'cst' structure defining 
% the CT images and the structure set. Make sure the matRad root directory with all its
% SUBDIRECTORIES is added to the Matlab search path.
clc,clear,close all
load('PROSTATE.mat');

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
% index vector depicting voxels in the CT cube that are covered by the corresponding VOI. In total, 17 structures are defined in the cst
cst

%%
% The fifth column represents meta parameters used for optimization such as the overlap priority, which can be specified in double presision. 
% A lower overlap priority indicates increased importance. In contrast, a higher overlap priority indicatets a strcture with lower importance. 
% The parameters alphaX and betaX depict the tissue's photon-radiosensitivity parameter of the linear quadratic model. These parameter
% are required for biological treatment planning using a variable RBE. Let's output the meta optimization parameter of the rectum, which is stored in the eighth row:
cst{8,2}
cst{8,5}

%%
% The sixth column contains optimization information such as objectives and constraints which are required to calculate the objective function value. 
% Please note, that multiple objectives/constraints can be defined for individual structures. As the rectum is an OAR, we have defined and
% squared overdosing objective so that it is considered to be expensive/costly for the optimizer delivering more than 50 Gy to the rectum. 
cst{8,6}

%% Treatment Plan
% The next step is to define your treatment plan labeld as 'pln'. This structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.

%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon. In this
% example we would like to use protons for treatment planning. Then, we need to define a treatment machine to correctly load the corresponding base data. Since we provide
% generic base data we set the machine to 'Genereric. By this means matRad will look for
% 'proton_Generic.mat' in our root directory and will use the data provided in there for dose calculation
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

%%
% Define the flavour of biological optimization for treatment planning along with the quantity that should be used for
% optimizaion. Possible values are (none: physical optimization; const_RBExD: constant RBE of 1.1; LEMIV_effect: 
% effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose. As we use protons, we follow here the clinical standard and use a constant 
% relative biological effectivness of 1.1. Therefore we set bioOptimization to const_RBExD
pln.bioOptimization = 'const_RBExD';   
                                       
%%
% Now we have to set some beam parameters. We can define multiple beam angles for the
% treatment and pass these to the plan as a vector. matRad will then interpret the vector as multiple beams.
% We define two opposing beams. For the first beam we set the gantry angle to 90 degree and the 
% corresponding couch angle to 0 degree. The second beam possess a gantry angle of 270 degree and a couch angle of 0 degree. 
% Furthermore, we want the lateral pencil beam spacing in x and y to be 3 mm in the iso-center plane which is perpendicular to the central beam ray.
% In total we are using 30 fractions. It is noteworthy that matRad is always optimizing the fraction dose.
pln.gantryAngles    = [90 270];
pln.couchAngles     = [0 0];
pln.bixelWidth      = 3;
pln.numOfFractions  = 30;

%%
% Obtain the number of beams and voxels from the existing variables and calculate the iso-center which is per default
% the mass of gravity of all target voxels.
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.voxelDimensions = ct.cubeDim;
pln.isoCenter       = ones(pln.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

%%
% Disable sequencing and direct aperture optimization, since we have a particle plan.
pln.runDAO        = 0;
pln.runSequencing = 0;

%%
% and ZACK! our treatment plan is ready. Lets have a look at it:
pln

%% Generatet Beam Geometry STF
% This acronym stands for steering file and comprises the complet beam geomtry along with 
% ray position, pencil beam positions and energies, source to axis distance (SAD) etc.
stf = matRad_generateStf(ct,cst,pln);

%%
% Let's display the beam geomtry information of the second beam
stf(2)

%% Dose Calculation
% Lets generate dosimetric information by pre-computing dose influence matrices for unit beamlet intensities. Having dose influences available allows then 
% later on inverse optimization. 
dij = matRad_calcParticleDose(ct,stf,pln,cst);

%% Inverse Optimizaiton for IMPT
% The goal of the fluence optimization is to find a set of bixel/spot weights which yield the best possible
% dose distribution according to the clinical objectives and constraints underlying the radiation treatment
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot the Resulting Dose Slice
% Let's plot the transversal iso-center dose slice
slice = round(pln.isoCenter(1,3)./ct.resolution.z);
figure
imagesc(resultGUI.RBExDose(:,:,slice)),colorbar,colormap(jet)

%%
% Now let's simulate a patient shift in y direction for both beams
stf(1).isoCenter(2) = stf(1).isoCenter(2) - 4;
stf(2).isoCenter(2) = stf(2).isoCenter(2) - 4;
pln.isoCenter       = reshape([stf.isoCenter],[3 pln.numOfBeams])';

%% Recalculate Plan
% Let's use the existing optimized pencil beam weights and recalculate the RBE weighted dose
resultGUI_isoShift = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);

%%  Visual Comparison of results
% Let's compare the new recalculation against the optimization result.
plane = 3;
doseWindow = [0 max([resultGUI.RBExDose(:); resultGUI_isoShift.RBExDose(:)])];

figure,title('original plan')
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.RBExDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);
figure,title('shifted plan')
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_isoShift.RBExDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);
 
absDiffCube = resultGUI.RBExDose-resultGUI_isoShift.RBExDose;
figure,title('absolute difference')
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],[],colorcube);

% Let's plot single profiles that are perpendicular to the beam direction
ixProfileY = round(pln.isoCenter(1,2)./ct.resolution.y);

profileOrginal = resultGUI.RBExDose(:,ixProfileY,slice);
profileShifted = resultGUI_isoShift.RBExDose(:,ixProfileY,slice);

figure,plot(profileOrginal,'LineWidth',2),grid on,hold on, 
       plot(profileShifted,'LineWidth',2),legend({'original profile','shifted profile'}),
       xlabel('mm'),ylabel('Gy(RBE)'),title('profile plot')
       
%% Quantitative Comparison of results
% Compare the two dose cubes using a gamma-index analysis. The gamma index
% is a composite quality distribution equaly taking into account a dose difference and a distance to
% criteria in order to reveal differences betweent two dose cubes.
% A gamma-index value of smaller than 1 indicates a successfull test and a
% value greater than 1 illustrates a failed test.

doseDifference   = 2;
distToAgreement  = 2;
n        = 1;

[gammaCube,gammaPassRateCell] = matRad_gammaIndex(...
    resultGUI_isoShift.RBExDose,resultGUI.RBExDose,...
    [ct.resolution.x, ct.resolution.y, ct.resolution.z],...
    slice,[doseDifference distToAgreement],n,'global',cst);


% Let's plot the gamma index histogram
figure,histogram(gammaCube(gammaCube>0),100),title('gamma index histogram')





