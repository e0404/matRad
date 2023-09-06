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

%% In this example we will show 
% (i) how to load patient data into matRad
% (ii) how to setup a photon dose calculation and 
% (iii) how to inversely optimize beamlet intensities
% (iv) how to visually and quantitatively evaluate the result
%global timed;
%timed = [];
%% Patient Data Import
% Let's begin with a clear Matlab environment. Then, import the TG119
% phantom into your workspace. The phantom is comprised of a 'ct' and 'cst'
% structure defining the CT images and the structure set. Make sure the 
% matRad root directory with all its subdirectories is added to the Matlab 
% search path.

matRad_rc; %If this throws an error, run it from the parent directory first to set the paths
%%
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultAccChangeTol = 1e-6;
%%

load('PROSTATE.mat');

%%

%%

pln.radiationMode = 'protons';  
%pln.radiationMode = 'photons';
pln.machine       = 'Generic';

quantityOpt    = 'RBExD';                                     
modelName      = 'constRBE';  

%pln.propDoseCalc.calcLET = 0;

pln.numOfFractions         = 30;
%pln.propStf.gantryAngles   = [0:52:359];
pln.propStf.gantryAngles   = [90 270];
pln.propStf.couchAngles    = zeros(1,numel(pln.propStf.gantryAngles));
pln.propStf.bixelWidth     = 5;


pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

%%
% Enable sequencing and disable direct aperture optimization (DAO) for now.
% A DAO optimization is shown in a seperate example.
pln.propOpt.runSequencing = 1;
pln.propOpt.runDAO        = 0;

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen');

%%
% and et voila our treatment plan structure is ready. Lets have a look:
display(pln);


%% Generate Beam Geometry STF
% The steering file struct comprises the complete beam geometry along with 
% ray position, pencil beam positions and energies, source to axis distance (SAD) etc.
stf = matRad_generateStf(ct,cst,pln);

%%
% Let's display the beam geometry information of the 6th beam
%display(stf(6));

%% Dose Calculation
% Let's generate dosimetric information by pre-computing dose influence 
% matrices for unit beamlet intensities. Having dose influences available 
% allows subsequent inverse optimization.
%dij = matRad_calcPhotonDose(ct,stf,pln,cst);
dij = matRad_calcParticleDose(ct,stf,pln,cst);
%% Inverse Optimization for IMRT
% The goal of the fluence optimization is to find a set of beamlet/pencil 
% beam weights which yield the best possible dose distribution according to
% the clinical objectives and constraints underlying the radiation 
% treatment. Once the optimization has finished, trigger once the GUI to 
% visualize the optimized dose cubes.
%% Paretooptimization
% The goal of this step is to define a grid of penalty values that
% are then evaluated using the matRad_paretoGeneration method
% The VOI and their respective penalties are defined in the following way
% It is possible to have more than one objective function per VOI
% penVal stores the Grid which is then passed on. penGrid contains an
% version easier to visualize, however harder to loop over
%% Assign objectives
cst{1,6}{1} = struct(DoseObjectives.matRad_EUDMin(100,0,8));
cst{8,6}{1} = struct(DoseObjectives.matRad_EUDMin(100,0,3));
cst{9,6}{1} = struct(DoseObjectives.matRad_MeanDose(100,0));
cst{6,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(61.2,78.2));
cst{7,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(50.4,64.4));
cst{1,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(0,58));
cst{9,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(0,58));
cst{8,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(0,58));
%%
cst{5,5}.Priority = 99;
%% Pareto optimization - comment out to run
%{
tic;            
retStruct = matRad_ParetoOptimization(dij,cst,pln,100);
toc
%}
%load('DataFull100Points.mat')
retStruct.weights = retStruct.weights*30
%%
%%
figure
plot(retStruct.errors,linewidth = 2)
xlabel('Iteration')
ylabel('Error [a.u]')
grid on
%% calculate example plans and plot their DVH and plan slices

slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
%% DVH calculation
resultGUIDoseGrid = matRad_calcCubesDoseGrid(retStruct.weights(:,1),dij)
resultGUIDoseGrid2 = matRad_calcCubesDoseGrid(retStruct.weights(:,2),dij)
resultGUIDoseGrid3 = matRad_calcCubesDoseGrid(retStruct.weights(:,5),dij)


%%
cstDVH = retStruct.modcst;
cstDVH(10,:) = [];
cstDVH(5,:) = [];
cstDVH(4,:) = [];
cstDVH(3,:) = [];
cstDVH(2,:) = [];
%%
dvh = matRad_calcDVH(cstDVH,resultGUIDoseGrid.RBExD);
dvh2 = matRad_calcDVH(cstDVH,resultGUIDoseGrid2.RBExD);
dvh3 = matRad_calcDVH(cstDVH,resultGUIDoseGrid3.RBExD);

figure
matRad_showDVH(dvh,cstDVH,pln);
matRad_showDVH(dvh2,cstDVH,pln,2);
matRad_showDVH(dvh3,cstDVH,pln,3);
%%
matlab2tikz('width', '\fwidth', 'height', '\fheight')
%% Plot slices
resultGUI = matRad_calcCubes(retStruct.weights(:,1),dij)
resultGUI2 = matRad_calcCubes(retStruct.weights(:,2),dij)
resultGUI3 = matRad_calcCubes(retStruct.weights(:,5),dij)
%%
cstI  = cst;
cstI(10,:) = [];
cstI(5,:) = [];
cstI(4,:) = [];
cstI(3,:) = [];
cstI(2,:) = [];

%%
figure
[hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(gca(),ct,cstI,1,resultGUI.RBExD,3,slice)
for j= 1: numel(hContour)
    hContour{j}(1).LineWidth = 2;
end
zoom(1.4)
set(gca(),'xlabel',[])
set(gca(),'ylabel',[])
set(gca(),'title',[])
set(gca(),'xtick',[])
set(gca(),'ytick',[])

%{
figure
[hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(gca(),ct,cstI,1,resultGUI2.RBExD,3,slice)
for j= 1: numel(hContour)
    hContour{j}(1).LineWidth = 2;
end
zoom(1.4)

figure
[hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(gca(),ct,cstI,1,resultGUI3.RBExD,3,slice)
for j= 1: numel(hContour)
    hContour{j}(1).LineWidth = 2;
end
zoom(1.4)
%}

matlab2tikz('width', '\fwidth', 'height', '\fheight')

%% For plot obtained from navigation 
figure
[hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(gca(),ct,cstI,1,resultGUI1.physicalDose,3,slice)
for j= 1: numel(hContour)
    hContour{j}(1).LineWidth = 2;
end
zoom(1.4)
%%
figure
plot(retStruct.errors)
%%
resultGUIDoseGridNav = matRad_calcCubesDoseGrid(resultGUI.w,dij)
dvhNav = matRad_calcDVH(cstDVH,resultGUIDoseGridNav.RBExD);
figure
matRad_showDVH(dvhNav,cstDVH,pln);
%%
matRad_plotParetoSurface(retStruct)
%%
matRad_plotPenaltyGrid(retStruct.penGrid)

%%
matRad_UIInterpolation(retStruct,dij,pln,ct,matRad_setOverlapPriorities(cst),retStruct.optiProb)