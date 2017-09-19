%% Example: Carbon Ion Treatment Plan
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
% Let's check the two structures, we have just imported. 'ct' comprises the ct cube 
% along with some meta information describing the ct cube (cube dimensions, resolution, number of CTs (usefull if
% we want to plan on a 4D ct).
ct

%%
% The 'cst' cell array defines volumes of interest/structures. Each row belongs to one certain VOI, whereas each column holds
% a different information. In particular, the second and third column show the name and the type of the structure. The fourth column depicts a linear 
% index list refering to the ct.cube and defines the actual membership. In total, 17 structures are defined in the cst
cst
%%
% The fifth column represents meta parameters such as the overlap priority. A lower overlap priority indicates increased importance. In contrast,
% a higher overlap priority indicatets a strcture with lower importance. The parameters alphaX and betaX depict the tissue's photon-radiosensitivity parameter
% which are required for biological treatment planning. Since the PTV, is stored in the tenth column, lets output the PTTVs meta information
cst{10,5}

%%
% The sixth column contains optimization information such as objective and constraints. Please note, that multiple objectives/constraints can be defined for a 
% individual structure
cst{10,6}

%% Now lets increase the importance of the square deviation objective to 500
cst{10,6}.penalty = 500;
%% Treatment Plan
% The next step is to define your treatment plan labeld as 'pln'. This structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.
%%
% First of all, we need to define what kind of radiation modality we would like to use.
% Then, we need to define a treattment machine the base data is associated with. Since we provide
% generic base data we set the machine to 'Genereric. By this means the the sofware will look for
% 'proton_Generic.mat' in our root directory and will use the data provided in there for dose calculation
pln.radiationMode = 'carbon';
pln.machine       = 'Generic';

%%
% Define the biological base data that should be used for treatment planning along with the quantity that should be used for
% optimizaion. LEMIV_RBExD results in using tabulated biological data based on LEMIV and use the RBE-weighted dose for optimization
pln.bioOptimization = 'LEMIV_RBExD';

%%
% Now we have to further specify the treatment plan. We define one beam with a gantry angle of 300 degree and set the 
% corresponding couch angle to 0 degree. Furthermore, we want the lateral spot spacing in x and y to be 3mm in the iso-center.
% In total we are using 30 fractions. It is noteworthy that matRad is always optimizing the fraction dose.
pln.gantryAngles    = 300;
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
% Now, let's have a final look at our treatment plan.
pln

%% Generatet the Beam Geometry
% This acronym stands for Steering File and comprises the beam geomtry along with 
% the ray and pencil beam positions
stf = matRad_generateStf(ct,cst,pln);

%% 
stf

%% Dose Calculation
% Calculate dose influence matrix for unit pencil beam intensities. 
dij = matRad_calcParticleDose(ct,stf,pln,cst);

%% Inverse Optimizaiton for IMPT
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot the Resulting Dose Slice
% Let's plot the transversal iso-center dose slice
slice = round(pln.isoCenter(3)./ct.resolution.z);
figure,
imagesc(resultGUI.RBExDose (:,:,slice)),colorbar

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
figure,
subplot(121)
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.RBExDose,plane,slice,[],[],colorcube);
subplot(122)
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_tissue.RBExDose,plane,slice,[],[],colorcube);

%% 
% At this point we would like to see the absolute difference of the optimization and the 
% recalculation. 
absDiffCube = resultGUI.RBExDose-resultGUI_tissue.RBExDose;
figure,[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],[],colorcube);




