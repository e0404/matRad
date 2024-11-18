%% Example: Proton Treatment Plan with Manipulated CT values including fine
%           sampling algorithm
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In this example we will show 
% (i) how to load patient data into matRad
% (ii) how to setup a proton dose calculation 
% (iii) how to inversely optimize the pencil beam intensities directly from command window in MATLAB.
% (iv) how to re-optimize a treatment plan
% (v) how to manipulate the CT cube by adding noise to the cube 
% (vi) how to recalculate the dose considering the manipulated CT cube and the previously optimized pencil beam intensities
% (vii) how to compare the two results

%% set matRad runtime configuration
matRad_rc; %If this throws an error, run it from the parent directory first to set the paths

%% Patient Data Import
% Let's begin with a clear Matlab environment and import the prostate 
% patient into your workspace.
load('PROSTATE.mat');

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.
pln.radiationMode                   = 'protons';           
pln.machine                         = 'generic_MCsquare'; %Use the base data fitted to MC here
pln.numOfFractions                  = 30;  
pln.bioModel                        = 'constRBE';
pln.multScen                        = 'nomScen';

%Geometry Settings
pln.propStf.gantryAngles            = [90 270];
pln.propStf.couchAngles             = [0 0];
pln.propStf.bixelWidth              = 5;
pln.propStf.longitudinalSpotSpacing = 5;
pln.propStf.numOfBeams              = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter               = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Optimization Settings
pln.propOpt.quantityOpt             = 'RBExDose'; 

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
%We do a Monte Carlo Dose calculation here to demonstrate how long an MC
%simulation on pencil-beam basis will take. If you just want to get through
%the example, feel free to use analytical dose calculation instead by
%uncommenting the first line or explicitly switch engine to HongPB

pln.propDoseCalc.engine = 'MCsquare';
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% Inverse Optimization for IMPT
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Calculate quality indicators 
resultGUI = matRad_planAnalysis(resultGUI,ct,cst,stf,pln);
ixRectum       = 1;
disp(resultGUI.qi(ixRectum).D_5);

%%
% Let's change the optimization parameter of the rectum in such a way that it
% will be better spared. We increase the penalty and lower the threshold 
% of the squared overdose objective function. Afterwards we re-optimize 
% the treatment plan and evaluate dose statistics one more time.

objective = cst{ixRectum,6}{1}; %This gives a struct
objective = matRad_DoseOptimizationFunction.createInstanceFromStruct(objective); %Now we turn it into a class
objective = objective.setDoseParameters(40); %We can simply call this function to change the/all dose parameter(s)
cst{ixRectum,6}{1} = struct(objective); % We put it back as struct for storage & compatability

cst{ixRectum,6}{1}.parameters{1} = 40;  % Change the reference dose
cst{ixRectum,6}{1}.penalty = 500; % Change the penalty
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
resultGUI = matRad_planAnalysis(resultGUI,ct,cst,stf,pln);

disp(resultGUI.qi(ixRectum).D_5);

%% Plot the Resulting Dose Slice
% Let's plot the transversal iso-center dose slice
slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
slice = slice(3);
figure
imagesc(resultGUI.RBExDose(:,:,slice)),colorbar, colormap(jet)

%% Add Range Uncertainty
% Now let's manually simulate a range undershoot by scaling the relative 
% stopping power cube by 3.5% percent. For that to happen, we need to tell 
% matRad that it should not convert the HU-cube (ct.cubeHU) to RSP cube 
% implicitly, but directly use the RSP cube we provide in ct.cube
% Note that such uncertainty scenarios can also be computed by using the
% functionalities of matRad_multScen

pln.propDoseCalc.useGivenEqDensityCube = true;
ct_manip         = ct;
ct_manip.cube{1} = 1.035*ct_manip.cube{1};

%% Recalculate Plan with MC square
% Let's use the existing optimized pencil beam weights and recalculate the RBE weighted dose
resultGUI_noise = matRad_calcDoseForward(ct_manip,cst,stf,pln,resultGUI.w);

%% Recalculate Plan with analytical fine sampling algorithm
% Again use the existing optimized pencil beam weights and recalculate the 
% RBE weighted dose, now using the fine sampling algorithm instead of MC

% pln.propDoseCalc.fineSampling stores parameters defining the fine 
% sampling simulation

pln.propDoseCalc.engine = 'SubsamplingPB';

pln.propDoseCalc.fineSampling.method = 'fitCircle';
pln.propDoseCalc.fineSampling.sigmaSub = 2; %mm
pln.propDoseCalc.fineSampling.N = 2;
% method for weight calculation, availabe methods:
%   'russo'
%   'fitCircle', supports N = 2,3 and 8
%   'fitSquare', supports N = 2 and 3

% pln.propDoseCalc.fineSampling.N = n sets the number of used fine
% sampling sub beams, default is N = 21
% parameter to modify number of calculated FS sub beams
%   'russo',        total number of beams = N^2
%   'fitCircle',    total number of beams = (2*N + 1)^2
%   'fitSquare',    total number of beams = (2^N - 1) * 6 + 1

% pln.propDoseCalc.fineSampling.sigmaSub = s set the Gaussian standard
% deviation of the sub Gaussian beams, only used when fine sampling
% method 'russo' is selected', default is s = 1;

% Direct call for fine sampling dose calculation:
resultGUI_FS = matRad_calcDoseForward(ct,cst,stf,pln,resultGUI.w);

%%  Visual Comparison of results using the "compareDose" helper function
matRad_compareDose(resultGUI_noise.RBExDose,resultGUI.RBExDose,ct,cst);
matRad_compareDose(resultGUI_FS.RBExDose,resultGUI.RBExDose,ct,cst);
