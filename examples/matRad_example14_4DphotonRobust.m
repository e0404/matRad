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

%% Import 4D CT
metadata.resolution = [3 3 3]; % [mm]
[ct,cst] = matRad_importMultipleDicomCt('/examples/dicom/patient_2',metadata);
clear 'metadata';

%% Instantiate elastic registration
metadata.nItera = 100;
metadata.dvfType = 'pull';
register = matRad_ElasticImageRegistration(ct,cst,1,metadata);
clear 'metadata';

%% Calculate deformation vector field
[ct,cst] = register.calcDVF();

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

%% Propagate contourns
%[ct,cst] = register.propContours();

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

%Activate 4D Optimization
%pln.propOpt.scen4D = 'all';
resultGUInominal = matRad_fluenceOptimization(dij,cst,pln);

% add resultGUInominal dose cubes to resultGUI structure to allow the visualization in the GUI
resultGUI = resultGUInominal;

%% Indicator calculation and show DVH and QI
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUInominal);

%% Trigger robust optimization
% Make the objective to a composite worst case objective
cst{ixCTV,6}{1}.robustness  = 'COWC';
cst{ixCTV,6}{2}.robustness  = 'COWC';
%cst{ixOAR,6}{1}.robustness  = 'COWC';

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

%% Indicator calculation and show DVH and QI
[dvh_robust,dqi_robust] = matRad_indicatorWrapper(cst,pln,resultGUIrobust);

%% calc 4D dose
totalPhaseMatrix = ones(dij.totalNumOfBixels,ct.numOfCtScen)/ct.numOfCtScen;  % the total phase matrix determines a mapping what fluence will be delivered in the which phase
totalPhaseMatrix = bsxfun(@times,totalPhaseMatrix,resultGUIrobust.w);         % equally distribute the fluence over all fluences

[resultGUIrobust4D, timeSequence] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUIrobust,totalPhaseMatrix);

% add resultGUIrobust4D dose cubes to the existing resultGUI structure to allow the visualization in the GUI
resultGUI = matRad_appendResultGUI(resultGUI,resultGUIrobust4D,0,'robust');

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
