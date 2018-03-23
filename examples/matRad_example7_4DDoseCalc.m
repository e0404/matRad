%% 4D dose calculation workflow
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Here we wil show the workflow for 4D dose calculation
%  
% Prerequisite for 4D dose calculation: 
%
% (1) 4D-CT incl. vector fields  
%
% Attention: first CT is reference CT used for optimization and as
% reference phase for dose accumulation
%
% (2) External program to calculate time delivery information (makeLmdout)
% 
% Assumptions:
% Regular breathing is assumed, breathing period and start point
% adjustable
%
%
% This example includes also scripts for
%
% * many 4D dose calculations with random breathing periods, start of
% delivery, delivery files
% * variable RBE for protons
% * 3d WC optimization and dose recalculation considering motion
%
%  Some of these tools include patient specific information (e.g. target
%  and OAR names and corresponding number in cst), that need to be adjustet
%  accordingly
%
% Silke Ulrich 2018
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%% Treatment planning
% First we plan the treatment (alternatively an existing treatment plan can
% be imported)

clc,clear,close all

addpath('D:\Matrad\')
addpath('D:\Matrad\4Ddose')
addpath('D:\Matrad\tools')

load('Liver_DS221.mat')

%%
% meta information for treatment plan
pln.numOfFractions  = 30;
pln.radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln.machine         = 'HITfixedBL';

% beam geometry settings
pln.propStf.bixelWidth      = 3; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longSpotSpacing = 3;      % only relevant for HIT machine, not generic
pln.propStf.gantryAngles    = [210 320]; 
pln.propStf.couchAngles     = [0 0]; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

%optimization settings
pln.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

quantityOpt  = 'RBExD';     % options: physicalDose, effect, RBExD
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions

scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'       

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,scenGenType);

% generate steering file
stf = matRad_generateStf(ct,cst,pln);

param.subIx = cst{4,4}{1};

% dose calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst,param);

% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln);


%% post processing
%This step is necessary to remove beam spots with too few particles that
%could not be delivered
resultGUI = matRad_postprocessing(resultGUI, dij, pln, cst, stf) ; 

%% export Plan
% Number of beams files called PBP_0X_Plan01.xml are generated 
% the order of spot positions defines the delivery path - you can choose
% between
% 'stfMode'- order of spots as in stf file - line wise
% 'backforth' - first row is delivered from left to right, next right to
% left and so on
% 'TSP' (attention, very slow) the shortest path between all spots in each
% energy slice is calculated
%
plnExportFilename = 'Plan01';
matRad_export_HITXMLPlan_modified(plnExportFilename, pln, stf, resultGUI, 'backforth')  


%% makeLmdout
% now we need to call the external program makeLmdout (installed on sievert22) to calculate the time
% structure of the delivery - this is not possible directly in matlab
% call: "./makeLmdout -p PBP_0X_Plan01.xml -o D_0X_Plan01 -y v2015" for all
% beams X - the binary is available on radf1 and runs on Sievert22 at
% /home/bangertm/lmdout. also note the necessary changes to
% /home/bangertm/.bash_login in order to run the binary
% the name of the output file should be the same as for the plan file (for matRad_calc4dDose) 
% if you call "./makeLmdout -p PBP_0X_Plan01.xml -o D_0X_Plan01 -y v2015 -
% x RANDOMINTEGER" and you supply a random integer after the argument
% identifer -x you can set a different random seed for the time sequence
% generation and thereby simulate uncertainty in the delivery sequence
% afterwards copy both the plan file (PBP...) and the output file (D...) to the 4dDose
% folder in matrad

%% calc 4D dose
% make sure that the correct pln, dij and stf are loeaded in the workspace
[resultGUI, delivery] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI, plnExportFilename); 

% Plot the result in comparison to the static dose
slice = round(pln.isoCenter(1,3)./ct.resolution.z); 
figure 
subplot(2,2,1)
imagesc(resultGUI.RBExD(:,:,slice)),colorbar, colormap(jet); 
title('static dose distribution')
subplot(2,2,2)
imagesc(resultGUI.accRBExD(:,:,slice)),colorbar, colormap(jet); 
title('4D dose distribution')
subplot(2,2,3)
imagesc(resultGUI.RBExD(:,:,slice) - resultGUI.accRBExD(:,:,slice)) ,colorbar, colormap(jet); 
title('Difference')


%% Appendix: Available scripts for 4D dose calculation
% some useful scripts related to 4D dose calculation

%% calc4DDose_loop
% a loop to calculate 4D dose for different lmdout files, different motion
% periods or different start of dose delivery (called motionoffset).

FileName = plnExportFilename;
GTVName = 'pseudoGTV';
OARName = 'Leber-GTV';
count = 2; % number of Lmdout files; the count has to be appended to the 
           % file name following an '_' (compare line 44, matRad_readLmdout)
MOTION = 'linear'; % spends same time in each motion phase
phaseTimeDev = 0; % only if sampling of different times in each phase is used
accDose = calc4dDose_loop(ct, pln, dij, stf, cst, resultGUI, FileName, GTVName, OARName, count, MOTION, phaseTimeDev);

% also a simulation of motion variation (not assuming perfectly regular motion) is possible
