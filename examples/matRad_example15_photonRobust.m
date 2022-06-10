%% Example: 4D robust Treatment Planning with photons
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team.
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
% (i)   import a 4D CT into a multiscenario ct and cst struct
% (ii)  create a photon treatment plan with seven beams
% (iii) perform dose calculation on each 4D CT
% (iv)  perform first a fluence optimization on the first CT scenario and then secondly
%       another fluence optimization using the composite worst case approach
%       considering all 4D CTs
% (v)   visualise all individual dose scenarios
% (vi)  sample discrete scenarios from Gaussian uncertainty assumptions

%% Clear workspace and command line
clear
clc
close 'all'

%% Set matRad runtime configuration
matRad_rc
param.logLevel=1;

%% Import 3D CT
load('patient3_scen_1.mat');

%% Plot CT slice
if param.logLevel == 1
    
    figure('Renderer', 'painters', 'Position', [10 10 300*ct.numOfCtScen 400]);
    
    isocenter = matRad_getIsoCenter(cst,ct,0);
    
    for scen_iterator = 1:ct.numOfCtScen
        plane      = 1;
        slice      = round(isocenter(2)./ct.resolution.y);
        subplot(2,ct.numOfCtScen,scen_iterator); camroll(90);
        matRad_plotSliceWrapper(gca,ct,cst,scen_iterator,[],plane,slice);
        
        plane      = 3;
        slice      = round(isocenter(3)./ct.resolution.z);
        subplot(2,ct.numOfCtScen,scen_iterator+ct.numOfCtScen);
        matRad_plotSliceWrapper(gca,ct,cst,scen_iterator,[],plane,slice);
        
    end
    
end
clear  scen_iterator plane slice ans;

if (ct.numOfCtScen>1)
    f          = figure; title('individual scenarios'); camroll(90);
    plane      = 1;
    slice      = round(isocenter(2)./ct.resolution.y);
    numScen    = 1;
    matRad_plotSliceWrapper(gca,ct,cst,numScen,[],plane,slice);
    b             = uicontrol('Parent',f,'Style','slider','Position',[50,5,419,23],...
        'value',numScen, 'min',1, 'max',ct.numOfCtScen,'SliderStep', [1/(ct.numOfCtScen-1) , 1/(ct.numOfCtScen-1)]);
    b.Callback    = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,round(es.Value),[],plane,slice);
end
clear  numScen plane slice ans f b;

%% Create the VOI data for the phantom

% Body
ixBody = 1;
cst{ixBody,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{ixBody,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,40));
cst{ixBody,6}{1}.robustness  = 'none';

% CTV
p=42.56;
ixCTV = 6;
cst{ixCTV,3}  = 'TARGET';
cst{ixCTV,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{ixCTV,6}{1} = struct(DoseObjectives.matRad_MinDVH(400,p,99));
cst{ixCTV,6}{1}.robustness  = 'none';
cst{ixCTV,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(1600,p));
cst{ixCTV,6}{2}.robustness  = 'none';
%cst{ixCTV,6}{3} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));

% Ipsilateral Lung
cst{2,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{2,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,30,16));
cst{2,6}{1}.robustness  = 'none';
cst{2,6}{2} = struct(DoseObjectives.matRad_MaxDVH(300,20,20));
cst{2,6}{2}.robustness  = 'none';
%cst{2,6}{3} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));

% Heart
cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(100,0));
cst{4,6}{1}.robustness  = 'none';

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
pln.numOfFractions         = 16;
pln.propStf.gantryAngles   = [11 42 73 104 135 309 340];
pln.propStf.couchAngles    = zeros(1,numel(pln.propStf.gantryAngles));
pln.propStf.bixelWidth     = 10;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 10; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 10; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 10; % [mm]

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);


%% save multi scenarios to plan

% retrieve the nominal scenario for dose calculation and optimziation
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
resultGUInominal = matRad_fluenceOptimization(dij,cst,pln);

% add resultGUInominal dose cubes to resultGUI structure to allow the visualization in the GUI
resultGUI = resultGUInominal;

%% Indicator calculation and show DVH and QI
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUInominal);

%% Define optimization scenarios
% retrieve 7 worst case scenarios for dose calculation and optimziation

multScen = matRad_multScen(ct,'wcScen');
multScen.wcFactor=1.0;
multScen.numOfShiftScen = [3 3 3];
multScen.shiftSD = [4 6 8];
multScen.includeNomScen=true;

%% copy nominal plan and save multi scenarios to plan
pln_robust=pln;
pln_robust.multScen = multScen;

%% Generate Beam Geometry STF
stf_robust = matRad_generateStf(ct,cst,pln_robust);

%% Dose Calculation
dij_robust = matRad_calcPhotonDose(ct,stf_robust,pln_robust,cst);

%% Trigger robust optimization
% Make the objective to a composite worst case objective
cst{ixCTV,6}{1}.robustness  = 'COWC';
cst{ixCTV,6}{2}.robustness  = 'COWC';

resultGUIrobustCOWC = matRad_fluenceOptimization(dij_robust,cst,pln_robust);

% add resultGUIrobust dose cubes to the existing resultGUI structure to allow the visualization in the GUI
resultGUI = matRad_appendResultGUI(resultGUI,resultGUIrobustCOWC,0,'COWC');

%% Setting robust plan for robustness analysis
resultGUIrobust=resultGUIrobustCOWC;

%% Indicator calculation and show DVH and QI
[dvh_robust,dqi_robust] = matRad_indicatorWrapper(cst,pln,resultGUIrobust);

%% Define sampling parameters
% select structures to include in sampling; leave empty to sample dose for all structures
% sampling does not know on which scenario sampling should be performed
structSel = {};

% Random scenarios
multScen = matRad_multScen(ct,'rndScen'); % 'impSamp' or 'wcSamp'
multScen.probDist = 'equalProb';
multScen.shiftSD = [4 6 8];
multScen.shiftGenType = 'sampled';
multScen.shiftCombType = 'combined';
multScen.numOfShiftScen = 25 * ones(3,1);
multScen.numOfRangeShiftScen = 25;
multScen.rangeRelSD=0;
multScen.rangeAbsSD=0;
multScen.scenCombType = 'combined';
multScen.includeNomScen=true;

%% Perform sampling for nominal optimization results
[caSamp, mSampDose, plnSamp, resultGUInomScen,resultGUIsampledScen] = matRad_sampling(ct,stf,cst,pln,resultGUI.w,structSel,multScen);

%% Perform sampling analysis
phaseProb = ones(1,ct.numOfCtScen)/ct.numOfCtScen;
varargin.robustnessCriterion = [10 10];
varargin.slice = round(isocenter(3)./ct.resolution.z);
[cstStat, resultGUISamp, meta] = matRad_samplingAnalysis(ct,cst,plnSamp,caSamp, mSampDose, resultGUInomScen,phaseProb,varargin);

%% Create an mean dose interactive plot to slide through axial slices
quantityMap='meanCubeW';
plane      = 3;
doseWindow = [0 max([max(resultGUISamp.(quantityMap)(:)) p/pln.numOfFractions*0.1])];
maxDose       = max([max(resultGUISamp.(quantityMap)(:)) p/pln.numOfFractions*0.1]);
doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
f = figure;
title([quantityMap 'for nominal optimization results']);
set(gcf,'position',[10,10,550,400]);
numSlices = ct.cubeDim(3);

slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(quantityMap),plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]');
b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
    'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(quantityMap),plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]');


%% Create an std dose interactive plot to slide through axial slices
quantityMap='stdCubeW';
plane      = 3;
doseWindow = [0 max([max(resultGUISamp.(quantityMap)(:)) p/pln.numOfFractions*0.1])];
maxDose       = max([max(resultGUISamp.(quantityMap)(:)) p/pln.numOfFractions*0.1]);
doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
f = figure;
title([quantityMap 'for nominal optimization results']);
set(gcf,'position',[10,10,550,400]);
numSlices = ct.cubeDim(3);

slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(quantityMap),plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]');
b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
    'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(quantityMap),plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]');

%% Multi-scenario dose volume histogram (DVH)
figure,set(gcf,'Color',[1 1 1],'position',[10,10,600,400]);
matRad_showDVHFromSampling(caSamp,cst,plnSamp,[1:plnSamp.multScen.totNumScen],[0 p/pln.numOfFractions*1.6],'trustband',1);
title('Multi-scenario DVH for nominal optimization results');

%% Perform sampling for robust optimization results
[caSampRob, mSampDoseRob, plnSampRob, resultGUIRobNomScen,resultGUIsampledScenRob] = matRad_sampling(ct,stf_robust,cst,pln_robust,resultGUIrobust.w,structSel,multScen);

%% Perform sampling analysis
phaseProb = ones(1,ct.numOfCtScen)/ct.numOfCtScen;
varargin.robustnessCriterion = [10 10];
varargin.slice = round(isocenter(3)./ct.resolution.z);
[cstStatRob, resultGUISampRob, metaRob] = matRad_samplingAnalysis(ct,cst,plnSampRob,caSampRob, mSampDoseRob, resultGUIRobNomScen,phaseProb,varargin);

%% Create an mean dose interactive plot to slide through axial slices
quantityMap='meanCubeW';
plane      = 3;
doseWindow = [0 max([max(resultGUISampRob.(quantityMap)(:)) p/pln.numOfFractions*0.1])];
maxDose       = max([max(resultGUISampRob.(quantityMap)(:)) p/pln.numOfFractions*0.1]);
doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
f = figure;
title([quantityMap 'for robust optimization results']);
set(gcf,'position',[10,10,550,400]);
numSlices = ct.cubeDim(3);

slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISampRob.(quantityMap),plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]');
b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
    'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISampRob.(quantityMap),plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]');

%% Create an std dose interactive plot to slide through axial slices
quantityMap='stdCubeW';
plane      = 3;
doseWindow = [0 max([max(resultGUISampRob.(quantityMap)(:)) p/pln.numOfFractions*0.1])];
maxDose       = max([max(resultGUISampRob.(quantityMap)(:)) p/pln.numOfFractions*0.1]);
doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
f = figure;
title([quantityMap 'for robust optimization results']);
set(gcf,'position',[10,10,550,400]);
numSlices = ct.cubeDim(3);

slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISampRob.(quantityMap),plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]');
b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
    'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISampRob.(quantityMap),plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]');

%% Multi-scenario dose volume histogram (DVH)
figure,set(gcf,'Color',[1 1 1],'position',[10,10,600,400]);
matRad_showDVHFromSampling(caSampRob,cst,plnSampRob,[1:plnSampRob.multScen.totNumScen],[0 p/pln.numOfFractions*1.6],'trustband',1);
title('Multi-scenario DVH for robust optimization results');

%% Perform price of robustness analysis in nominal scenario
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
[resultGUISampRob] = matRad_priceOfRobustnessIndex(resultGUISampRob,resultGUInominal.(quantityOpt),resultGUIRobNomScen.(quantityOpt),ct,cst,[],[],[],[-5 5],'relative',slice);

%% Perform price of robustness analysis using mean dose
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
[resultGUISampRob] = matRad_priceOfRobustnessIndex(resultGUISampRob,resultGUISamp.meanCubeW,resultGUISampRob.meanCubeW,ct,cst,[],[],[],[-5 5],'relative',slice);
