%% 4D dose calculation workflow
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Here we wil show the workflow for 4D dose calculation
%  
% You can find the data (LiverDS221 mat data and treatment plan PBP... and lmdout output file (D...)) used for this example here:  \\radfsjulia\E040\E0404\DataSilke_4D\Liver007
%
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
% Silke Ulrich 2018
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%% Treatment planning
% First we plan the treatment (alternatively an existent treatment plan can
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
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longSpotSpacing = 5;      % only relevant for HIT machine, not generic
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

% dose calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst,false);

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
matRad_export_HITXMLPlan_modified('Plan01_5new',  pln, stf, resultGUI, 'backforth')  

%% makeLmdout
% now we need to call the external program makeLmdout (installed on sievert22) to calculate the time
% structure of the delivery - this is not possible directly in matlab
% Aufruf: makeLmdout -p PBP_0X_Plan01.xml -o D_0X_Plan01 -y v2015 für jeden
% beam X
% the name of the output file should be the same as for the plan file (for matRad_calc4dDose) 
% copy both the plan file (PBP...) and the output file (D...) to the 4dDose
% folder in matrad

%% calc 4D dose
% make sure that the correct pln, dij and stf are loeaded in the workspace
[resultGUI, delivery] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI, 'Plan01_5new'); 

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


