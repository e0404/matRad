%% Example: 4D robust Treatment Planning with photons
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In this example we will  
% (i)   create a small artifical phantom
% (ii)  add a movement to the phantom to create a 4D CT
% (iii) create a photon treatment plan with two beams
% (iv)  perform dose calculation on each 4D CT
% (iv)  perform first a fluence optimization on the first CT scenario and then secondly
%       another fluence optimization using the composite worst case paradigm 
%       considering all 4D CTs
% (v)   visualise all individual dose scenarios 
% (vi)  sample discrete scenarios from Gaussian uncertainty assumptions

%% set matRad runtime configuration
matRad_rc

%% Create an artifiical CT image series
xDim = 150;
yDim = 150;
zDim = 50;

ct.cubeDim      = [xDim yDim zDim];
ct.resolution.x = 2; % mm
ct.resolution.y = 2; % mm
ct.resolution.z = 3; % mm
ct.numOfCtScen  = 1;
 
% create an ct image series with zeros - it will be filled later
ct.cubeHU{1} = ones(ct.cubeDim) * -1024;

%% Create the VOI data for the phantom
% Now we define two structures for the phantom 
ixOAR = 1;
ixPTV = 2;

% define general VOI properties
cst{ixOAR,1} = 0;
cst{ixOAR,2} = 'contour';
cst{ixOAR,3} = 'OAR';
cst{ixPTV,1} = 1;
cst{ixPTV,2} = 'target';
cst{ixPTV,3} = 'TARGET';
 
% define optimization parameter for both VOIs
cst{ixOAR,5}.TissueClass = 1;
cst{ixOAR,5}.alphaX      = 0.1000;
cst{ixOAR,5}.betaX       = 0.0500;
cst{ixOAR,5}.Priority    = 2;           % overlap priority for optimization - a higher number corresponds to a lower priority
cst{ixOAR,5}.Visible     = 1;

cst{ixOAR,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10,30));

cst{ixPTV,5}.TissueClass = 1;
cst{ixPTV,5}.alphaX      = 0.1000;
cst{ixPTV,5}.betaX       = 0.0500;
cst{ixPTV,5}.Priority    = 1;           % overlap priority for optimization - a lower number corresponds to a higher priority
cst{ixPTV,5}.Visible     = 1; 

cst{ixPTV,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(50,60));


%% Lets create a cubic phantom
% first define the dimensions of the OAR
cubeHelper = zeros(ct.cubeDim);
xLowOAR    = round(xDim/2 - xDim/6);
xHighOAR   = round(xDim/2 + xDim/6);
yLowOAR    = round(yDim/2 - yDim/6);
yHighOAR   = round(yDim/2 + yDim/6);
zLowOAR    = round(zDim/2 - zDim/4);
zHighOAR   = round(zDim/2 + zDim/4);

for x = xLowOAR:1:xHighOAR
   for y = yLowOAR:1:yHighOAR
      for z = zLowOAR:1:zHighOAR
         cubeHelper(x,y,z) = 1;
      end
   end
end
      
% extract the linear voxel indices and save it in the cst
cst{ixOAR,4}{1} = find(cubeHelper);

% second the PTV
cubeHelper = zeros(ct.cubeDim);
radiusPTV = xDim/14;
for x = 1:xDim
   for y = 1:yDim
      for z = 1:zDim
         currPost = [x y z] - round([ct.cubeDim./2]);
         if  sqrt(sum(currPost.^2)) < radiusPTV
            cubeHelper(x,y,z) = 1;
         end
      end
   end
end

% extract the linear voxel indices and save it in the cst
cst{ixPTV,4}{1} = find(cubeHelper);

% assign relative electron densities
vIxOAR = cst{ixOAR,4}{1};
vIxPTV = cst{ixPTV,4}{1};

ct.cubeHU{1}(vIxOAR) = 300; % assign HU of soft tissue
ct.cubeHU{1}(vIxPTV) = 0;   % assign HU of water

%% add motion to the phantom and artificially create a 4D CT with vector fields
amplitude    = [5 0 0]; % [voxels]
numOfCtScen  = 5;
motionPeriod = 2.5; % [s] 

[ct,cst] = matRad_addMovement(ct, cst,motionPeriod, numOfCtScen, amplitude,1);

% show the deformation vector field
slice = 25; % select a specific slice and to plot the vector field
[a,xDim,yDim,zDim] = size(ct.dvf{1});

[mX,mY]      = meshgrid(1:xDim,1:yDim);

figure,
for ctPhase = 1:ct.numOfCtScen 
   clf;
   xVectorField = squeeze(ct.dvf{ctPhase}(1,:,:,slice));  % retrieve the deformation vector field in x-direction of slice 25
   yVectorField = squeeze(ct.dvf{ctPhase}(2,:,:,slice));  % retrieve the deformation vector field in y-direction of slice 25
   quiver(mX,mY,yVectorField,xVectorField); title(['deformation vector field of phase ' num2str(ctPhase)]),
   set(gca,'XLim',[70 80]);set(gca,'YLim',[70 80]);
   % flip y axis to be consistent with the previous plot
   ax = gca; set(ax,'YDir','reverse');
   pause(0.5);  
end

% magnitude of the vector field should change over ct scenarios
% vector field refers to the first (initial) ct scenario
% investigate in the difference of 4D dose accumulatino


% clear helper variables to get clean workspace
clear x y z xDim yDim zDim xHighOAR xLowOAR xHighOAR yHighOAR yLowOAR zHighOAR zLowOAR vIxOAR vIxPTV cubeHelper currPost radiusPTV
%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines the most
% important cornerstones of your treatment plan.
%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon. In this
% example we would like to use protons for robust treatment planning. Next, we
% need to define a treatment machine to correctly load the corresponding 
% base data. matRad features generic base data in the file
% 'carbon_Generic.mat'; consequently the machine has to be set accordingly
pln.radiationMode = 'photons';            
pln.machine       = 'Generic';

%%
% Define the biological optimization model for treatment planning along
% with the quantity that should be used for optimization. Possible model values 
% are:
% 'none':     physical optimization;
% 'constRBE': constant RBE of 1.1; 
% 'MCN':      McNamara-variable RBE model for protons; 
% 'WED':      Wedenberg-variable RBE model for protons
% 'LEM':      Local Effect Model 
% and possible quantityOpt are 'physicalDose', 'effect' or 'RBExD'.
% As we use protons, we use a constant RBE of 1.1.
modelName    = 'none';
quantityOpt  = 'physicalDose';   

%%
% The remaining plan parameters are set like in the previous example files
pln.numOfFractions        = 20;
pln.propStf.gantryAngles  = [0 90];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% retrieve 9 worst case scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen');                                         

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
dij = matRad_calcPhotonDose(ct,stf,pln,cst);

%% Inverse Optimization  for IMPT based on RBE-weighted dose
% The goal of the fluence optimization is to find a set of bixel/spot 
% weights which yield the best possible dose distribution according to the
% clinical objectives and constraints underlying the radiation treatment.
% 

%Activate 4D Optimization
%pln.propOpt.scen4D = 'all';
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Trigger robust optimization
% Make the objective to a composite worst case objective
cst{ixPTV,6}{1}.robustness  = 'COWC';
cst{ixOAR,6}{1}.robustness  = 'COWC';

%Activate 4D Optimization
pln.propOpt.scen4D = 'all';

% voxel wise worst case 'VWWC' does not work for 4D robust optimization

% parameters for stochastic optimization
%cst{ixPTV,6}.robustness  = 'STOCH';
%cst{ixOAR,6}.robustness  = 'STOCH';
%pln.multScen.scenProb = (1/ct.numOfCtScen) * ones(ct.numOfCtScen,1);  % assign probabilities to 4D scenarios

% % parameters for objective wise worst case
%cst{ixPTV,6}.robustness  = 'OWC';
%cst{ixOAR,6}.robustness  = 'OWC';

resultGUIrobust = matRad_fluenceOptimization(dij,cst,pln);

% add resultGUIrobust dose cubes to the existing resultGUI structure to allow the visualization in the GUI
resultGUI = matRad_appendResultGUI(resultGUI,resultGUIrobust,0,'robust');

%% calc 4D dose
totalPhaseMatrix = ones(dij.totalNumOfBixels,ct.numOfCtScen)/ct.numOfCtScen;  % the total phase matrix determines a mapping what fluence will be delivered in the which phase
totalPhaseMatrix = bsxfun(@times,totalPhaseMatrix,resultGUIrobust.w);         % equally distribute the fluence over all fluences

[resultGUIrobust4D, timeSequence] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUIrobust,totalPhaseMatrix); 

%% Visualize results

plane         = 3;
slice         = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
maxDose       = max([max(resultGUI.([quantityOpt])(:,:,slice)) max(resultGUIrobust.([quantityOpt])(:,:,slice))])+1e-4;
doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
figure,
subplot(121),matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.([quantityOpt '_' 'beam1'])      ,plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);title('conventional plan - beam1')
subplot(122),matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.([quantityOpt])      ,plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);title('conventional plan')

figure
subplot(121),matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIrobust.([quantityOpt '_' 'beam1']),plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);title('robust plan - beam1')
subplot(122),matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIrobust.([quantityOpt]),plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);title('robust plan')


figure
subplot(131),matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIrobust4D.([quantityOpt]),plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);title('robust plan')
subplot(132),matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIrobust4D.('accPhysicalDose'),plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);title('robust plan dose accumulation')

% create an interactive plot to slide through individual scnearios
f       = figure; title('individual scenarios');
numScen = 1;
maxDose       = max(max(resultGUIrobust.([quantityOpt '_' num2str(round(numScen))])(:,:,slice)))+0.2;
doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIrobust.([quantityOpt '_' num2str(round(numScen))]),plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);

[env,envver] = matRad_getEnvironment();
if strcmp(env,'MATLAB') || str2double(envver(1)) >= 5
    b = uicontrol('Parent',f,'Style','slider','Position',[50,5,419,23],...
        'value',numScen, 'min',1, 'max',pln.multScen.totNumScen,'SliderStep', [1/(pln.multScen.totNumScen-1) , 1/(pln.multScen.totNumScen-1)]);
    set(b,'Callback',@(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,round(get(es,'Value')),resultGUIrobust.([quantityOpt '_' num2str(round(get(es,'Value')))]),plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels));
end

%% Indicator calculation and show DVH and QI
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUIrobust);

%% Perform sampling
% select structures to include in sampling; leave empty to sample dose for all structures

% sampling does not know on which scenario sampling should be performed
structSel = {}; % structSel = {'PTV','OAR1'};
[caSamp, mSampDose, plnSamp, resultGUInomScen]          = matRad_sampling(ct,stf,cst,pln,resultGUI.w,structSel);
[cstStat, resultGUISamp, meta]                         = matRad_samplingAnalysis(ct,cst,plnSamp,caSamp, mSampDose, resultGUInomScen);

[caSampRob, mSampDoseRob, plnSampRob, resultGUInomScen] = matRad_sampling(ct,stf,cst,pln,resultGUIrobust.w,structSel);
[cstStatRob, resultGUISampRob, metaRob]                = matRad_samplingAnalysis(ct,cst,plnSampRob,caSampRob, mSampDoseRob, resultGUInomScen);

figure,title('std dose cube based on sampling - conventional')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.stdCube,plane,slice,[],[],colorcube,[],[0 max(resultGUISamp.stdCube(:))]);

figure,title('std dose cube based on sampling - robust')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISampRob.stdCube,plane,slice,[],[],colorcube,[],[0 max(resultGUISampRob.stdCube(:))]);



