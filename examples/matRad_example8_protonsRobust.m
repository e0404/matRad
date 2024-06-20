%% Example: Robust Treatment Planning with Protons
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
% (ii)  create a scanned proton treatment plan considering a constant RBE of 1.1
% (iii) we will enable dose calculation on nine selected worst case scenarios
% (iv)  robustly optimize the pencil beam intensities on all 9 dose scenarios 
%       using the composite worst case paradigm 
% (v)   visualise all individual dose scenarios 
% (vi)  sample discrete scenarios from Gaussian uncertainty assumptions

%% set matRad runtime configuration
matRad_rc

%% Create a CT image series with the Phantom Builder

ctDim = [120,120,60]; % x,y,z dimensions
ctResolution = [3,3,3]; % x,y,z the same here!

builder = matRad_PhantomBuilder(ctDim,ctResolution,1);

% Now we define three structures for the phantom 
objective1 = struct(DoseObjectives.matRad_SquaredDeviation(800,45));
objective2 = struct(DoseObjectives.matRad_SquaredOverdosing(400,0));

builder.addSphericalTarget('target',ctDim(1)/13,'objectives',struct(DoseObjectives.matRad_SquaredDeviation(100,60)));
builder.addSphericalOAR('OAR',ctDim(1)/15,'offset',[-10 10 0],'HU',-100,'objectives',struct(DoseObjectives.matRad_SquaredOverdosing(10,40)));
builder.addBoxOAR('contour',ctDim./2,'HU',0,'objectives',struct(DoseObjectives.matRad_SquaredOverdosing(5,20)));

%Keep indices for assignment of robustness objectives later
ixTarget = 1;
ixOAR = 2;
ixNT = 3;

[ct,cst] = builder.getctcst();
  
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
pln.radiationMode = 'protons';            
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
modelName    = 'constRBE';
quantityOpt  = 'RBExD';   

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

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% retrieve 9 worst case scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'wcScen');                                         

pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

%% Inverse Optimization  for IMPT based on RBE-weighted dose
% The goal of the fluence optimization is to find a set of bixel/spot 
% weights which yield the best possible dose distribution according to the
% clinical objectives and constraints underlying the radiation treatment.

resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Trigger robust optimization
% Make the objective to a composite worst case objective

ROBUST_OPT = {'COWC'}; %{'STOCH','VWWC','VWWC_INV','COWC','OWC','PROB'};

for ixRob = 1:numel(ROBUST_OPT)
   cst{ixTarget,6}{1}.robustness  = ROBUST_OPT{1,ixRob};
   cst{ixOAR,6}{1}.robustness     = ROBUST_OPT{1,ixRob};
   cst{ixNT,6}{1}.robustness      = ROBUST_OPT{1,ixRob};
   
   % add a max constraint
   %cst{ixOAR,6}{1,2} = struct(DoseConstraints.matRad_MinMaxDose([0 20],'voxel'));
   %cst{ixOAR,6}{1,2}.robustness   = 'COWC';
   
   resultGUIrobust = matRad_fluenceOptimization(dij,cst,pln);
   
   % combine resultGUI structures
   resultGUI = matRad_appendResultGUI(resultGUI,resultGUIrobust,0,['robust' ROBUST_OPT{1,ixRob}]);
   
end

% matRadGUI

%% Visualize results
plane      = 3;
slice      = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);

figure,matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.RBExD_beam1      ,plane,slice,[],[],colorcube,[],[0 max(resultGUI.RBExD_beam1(:))],[]);title('conventional plan - beam1')
figure,matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIrobust.RBExD_beam1,plane,slice,[],[],colorcube,[],[0 max(resultGUIrobust.RBExD_beam1(:))],[]);title('robust plan - beam1')

% create an interactive plot to slide through individual scnearios
f = figure;title('individual scenarios');
numScen = 1;doseWindow = [0 3.5];
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIrobust.(['RBExD_' num2str(round(numScen))]),plane,slice,[],[],colorcube,[],doseWindow,[]);

[env,envver] = matRad_getEnvironment();
if strcmp(env,'MATLAB') || str2double(envver(1)) >= 5
    b = uicontrol('Parent',f,'Style','slider','Position',[50,5,419,23],...
        'value',numScen, 'min',1, 'max',pln.multScen.totNumScen,'SliderStep', [1/(pln.multScen.totNumScen-1) , 1/(pln.multScen.totNumScen-1)]);
    set(b,'Callback',@(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIrobust.(['RBExD_' num2str(round(get(es,'Value')))]),plane,slice,[],[],colorcube,[],doseWindow,[]));
end

%% Indicator calculation and show DVH and QI
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUIrobust);
%% Perform sampling
% select structures to include in sampling; leave empty to sample dose for all structures
structSel = {}; % structSel = {'PTV','OAR1'};
[caSamp, mSampDose, plnSamp, resultGUInomScen]          = matRad_sampling(ct,stf,cst,pln,resultGUI.w,structSel);
[cstStat, resultGUISamp, meta]                         = matRad_samplingAnalysis(ct,cst,plnSamp,caSamp, mSampDose, resultGUInomScen);

[caSampRob, mSampDoseRob, plnSampRob, resultGUInomScen] = matRad_sampling(ct,stf,cst,pln,resultGUIrobust.w,structSel);
[cstStatRob, resultGUISampRob, metaRob]                = matRad_samplingAnalysis(ct,cst,plnSampRob,caSampRob, mSampDoseRob, resultGUInomScen);

figure,title('std dose cube based on sampling - conventional')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.stdCube,plane,slice,[],[],colorcube,[],[0 max(resultGUISamp.stdCube(:))]);

figure,title('std dose cube based on sampling - robust')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISampRob.stdCube,plane,slice,[],[],colorcube,[],[0 max(resultGUISampRob.stdCube(:))]);

