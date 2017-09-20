%% Example Photon Treatment Plan
% In this example we will show how to load some patient data on matrad and
% how to setup a simulation directly from command window in MatLab.

%% Import
% Let's begin importing examples ct and cst files in our workspace and 
% adding the path of the main folder.
addpath('..\')
load('TG119.mat');

%%
% We can check if we have all the information we need, calling the two
% imported files. 'ct' is a very simple structure containing just: the ct
% cube; its dimension; its resolution; and the number of scenes (usefull if
% we want to plan on a 4D ct).
ct

%%
% In 'cst' file we have the data of the different areas inside the ct.
% Apart from index, name and classification, we have the indeces
% representing the volume of that specific area in a cell in the fouth
% column.
cst
%%
% Fifth column represents meta parameters or tissue information. and
% options for tissue visualization.
cst{1,5}

%%
% Sixth column contains optimization information.
cst{1,6}

%% Plan
% Next step is to build our treatment plan. This has different parameters
% that we will set together and we usually refer to it as 'pln'.
%%
%%
% First of all, we need to set the some basic parameters about the machine.
% We add the chosen particle for treatment then we set the substruct
% machine as 'Generic' so in this way the sofware will look for
% 'photon_Generic.mat' in our folder and will assign those as machine data
% for this simulation.
pln.radiationMode = 'photons';
pln.machine = 'Generic';


%%
% Now we set some beam parameters. We can chose multiple angles for the
% treatment and pass these to the plan as a vector and then MatRad will now
% what to do. In this case, we chose 5 beams with 5 mm bixel width.
pln.gantryAngles = 0:72:359;
pln.couchAngles = [0, 0, 0, 0, 0];
pln.numOfBeams = 5;
pln.bixelWidth = 5;

%%
% Some physical parameters still need to be set. We repeat the isocenter
% coordinates for every entering beam.
pln.isoCenter = repmat([83, 83, 64],pln.numOfBeams,1);
pln.voxelDimensions = ct.cubeDim;
pln.numOfVoxels = prod(pln.voxelDimensions);

%%
% Last, some extra feature..
pln.bioOptimization = 0;
pln.runDAO = 0;
pln.runSequencing = 0;

%%
% and our Plan is ready!
pln

%% Stf
% This acronym stands for Steering File. There is a very simple way to
% generate this, we just call the following function:
stf = matRad_generateStf(ct,cst,pln);

%%
% And we get, for example for the third beam:
stf(1,3)

%% Dose Calculation
% In order to run a calculation, we simply create weights data (that will
% be corrected by optimization and call:
w = ones(1, sum([stf.numOfBixelsPerRay]));
resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,w);

%%
% Now let's try to take one beam out and repeate the calculation.
pln.gantryAngles = 0:72:220;
pln.couchAngles = [0, 0, 0, 0];
pln.numOfBeams = 4;
%%
w = ones(1, sum([stf.numOfBixelsPerRay]));
resultGUI_4beam = matRad_calcDoseDirect(ct,stf,pln,cst,w);

%% Plot
% Now you surely want to see what the result of this computation is.
% The simpliest way to do this is to call
figure
imagesc(resultGUI.physicalDose(:,:,stf(1,1).isoCenter(3)))
%%
% Or you can use a MatRad function that Plots not only a slice of the Dose 
% but it does it on the appropriate ct slice.
%%
% Let's say we would like to plot the slice that corresponds to the
% isocenter on the axial plane. We set
plane = 3;
slice = stf(1,1).isoCenter(3);

%%
% and get
addpath('tools\')
figure
[hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],[],colorcube)
figure
[hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_4beam.physicalDose,plane,slice,[],[],colorcube)

%% Comparisons
% At this point one would like to compare the simulations we made with each
% other, to see the differences and maybe check the validity of a certain
% prediction. We should now use another MatRad usefull function that will
% give us a map of the so called gamma-index (method described on Low et
% al. 1998).
%%
% We set some parameters before.
% 'criteria' indicates the level of accuracy we look for, in this case, we
% set them at 3% of the Dose and 1mm of distance. 'n' is the number of
% times we want to interpolate the matrix before the comparison. We suggest
% to use values from 1 to 3, you can set it higher at you own risk.
criteria = [3, 1];
n = 1;

%%
% Here we can call the function, after including folder in the path
[gammaCube,gammaPassRateCell] = matRad_gammaIndex(...
    resultGUI_4beam.physicalDose,resultGUI.physicalDose,...
    [ct.resolution.x, ct.resolution.y, ct.resolution.z],...
    pln.isoCenter(3),criteria,n,'global',cst);








