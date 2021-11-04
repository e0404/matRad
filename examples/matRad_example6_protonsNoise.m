%% Example: Proton Treatment Plan with Manipulated CT values including fine
%           sampling algorithm
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
% (iii) how to inversely optimize the pencil beam intensities directly from command window in MATLAB.
% (iv) how to re-optimize a treatment plan
% (v) how to manipulate the CT cube by adding noise to the cube 
% (vi) how to recalculate the dose considering the manipulated CT cube and the previously optimized pencil beam intensities
% (vii) how to compare the two results

%% Patient Data Import
% Let's begin with a clear Matlab environment and import the prostate 
% patient into your workspace.

matRad_rc; %If this throws an error, run it from the parent directory first to set the paths

load('PROSTATE.mat');

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.

pln.radiationMode           = 'protons';           
pln.machine                 = 'generic_MCsquare';
pln.numOfFractions          = 30;
pln.propOpt.bioOptimization = 'const_RBExD';     
pln.propStf.gantryAngles    = [90 270];
pln.propStf.couchAngles     = [0 0];
pln.propStf.bixelWidth      = 3;
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO          = 0;
pln.propOpt.runSequencing   = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
dij = matRad_calcParticleDoseMC(ct,stf,pln,cst);

%% Inverse Optimization for IMPT
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Calculate quality indicators 
[dvh,qi]       = matRad_indicatorWrapper(cst,pln,resultGUI);
ixRectum       = 1;
display(qi(ixRectum).D_5);

%%
% Let's change the optimization parameter of the rectum in such a way that it
% will be better spared. We increase the penalty and lower the threshold 
% of the squared overdose objective function. Afterwards we re-optimize 
% the treatment plan and evaluate dose statistics one more time.

objective = cst{ixRectum,6}{1}; %This gives a struct
objective = matRad_DoseOptimizationFunction.createInstanceFromStruct(objective); %Now we turn it into a class
objective = objective.setDoseParameters(40); %We can simply call this function to change the/all dose parameter(s)
cst{ixRectum,6}{1} = struct(objective); % We put it back as struct

cst{ixRectum,6}{1}.parameters{1} = 40;
cst{ixRectum,6}{1}.penalty = 500;
resultGUI               = matRad_fluenceOptimization(dij,cst,pln);
[dvh2,qi2]              = matRad_indicatorWrapper(cst,pln,resultGUI);
display(qi2(ixRectum).D_5);

%% Plot the Resulting Dose Slice
% Let's plot the transversal iso-center dose slice
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
figure
imagesc(resultGUI.RBExDose(:,:,slice)),colorbar, colormap(jet)

%%
% Now let's simulate a range undershoot by scaling the relative stopping power cube by 3.5% percent
ct_manip         = ct;
ct_manip.cubeHU{1} = 1.035*ct_manip.cubeHU{1};

%% Recalculate Plan with MC square
% Let's use the existing optimized pencil beam weights and recalculate the RBE weighted dose
resultGUI_noise = matRad_calcDoseDirectMC(ct_manip,stf,pln,cst,resultGUI.w);

%% Recalculate Plan with analytical fine sampling algorithm
% Again use the existing optimized pencil beam weights and recalculate the 
% RBE weighted dose, now using the fine sampling algorithm instead of MC

% pln.propDoseCalc.fineSampling stores parameters defining the fine 
% sampling simulation

pln.propDoseCalc.fineSampling.method = 'russo'; 
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

% Indirect call for fine sampling dose calculation:    
% dijFS = matRad_calcParticleDose(ct,stf,pln,cst,false);
% resultGUI_FS = matRad_calcCubes(resultGUI.w,dijFS);

% Direct call for fine sampling dose calculation:
resultGUI_FS = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);

%%  Visual Comparison of results
% Let's compare the new recalculation against the optimization result.
plane      = 3;
doseWindow = [0 max([resultGUI.RBExDose(:); resultGUI_noise.RBExDose(:)])];

figure,title('original plan')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.RBExDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);
figure,title('manipulated plan MC')
matRad_plotSliceWrapper(gca,ct_manip,cst,1,resultGUI_noise.RBExDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);
figure,title('manipulated plan FS')
matRad_plotSliceWrapper(gca,ct_manip,cst,1,resultGUI_FS.RBExDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);

% Let's plot single profiles along the beam direction
ixProfileY = round(pln.propStf.isoCenter(1,1)./ct.resolution.x);

profileOrginal = resultGUI.RBExDose(:,ixProfileY,slice);
profileNoiseMC   = resultGUI_noise.RBExDose(:,ixProfileY,slice);
profileNoiseFS   = resultGUI_FS.RBExDose(:,ixProfileY,slice);

figure,plot(profileOrginal,'LineWidth',2),grid on,hold on, 
       plot(profileNoiseMC,'LineWidth',2),plot(profileNoiseFS,'LineWidth',2),
       legend({'original profile','noise profile MC','noise profile FS'}),
       xlabel('mm'),ylabel('Gy(RBE)'),title('profile plot')
       
%% Quantitative Comparison of results
% Compare the two dose cubes using a gamma-index analysis.

doseDifference   = 2;
distToAgreement  = 2;
n                = 1;

[gammaCubeMC,gammaPassRateCellMC] = matRad_gammaIndex(...
    resultGUI_noise.RBExDose,resultGUI.RBExDose,...
    [ct.resolution.x, ct.resolution.y, ct.resolution.z],...
    [doseDifference distToAgreement],slice,n,'global',cst);

[gammaCubeFS,gammaPassRateCellFS] = matRad_gammaIndex(...
    resultGUI_FS.RBExDose,resultGUI.RBExDose,...
    [ct.resolution.x, ct.resolution.y, ct.resolution.z],...
    [doseDifference distToAgreement],slice,n,'global',cst);

