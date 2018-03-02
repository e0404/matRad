% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

% load patient data, i.e. ct, voi, cst

%load HEAD_AND_NECK
%load TG119.mat
%load PROSTATE.mat
%load LIVER.mat
load BOXPHANTOM.mat

%% initial visualization and change objective function settings if desired
matRadGUI

% meta information for treatment plan
pln.numOfFractions  = 30;
pln.radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln.machine         = 'Generic';

% beam geometry settings
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longSpotSpacing = 3;      % only relevant for HIT machine, not generic
pln.propStf.gantryAngles    = [0:72:359]; % [?] ;
pln.propStf.couchAngles     = [0 0 0 0 0]; % [?] ; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

pln.propOpt.minNrParticles  = 500000;

pln.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions

scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,scenGenType);

%% initial visualization and change objective function settings if desired
matRadGUI

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);

%% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst,false);
    %dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
elseif any(strcmp(pln.radiationMode,{'protons','helium','carbon'}))
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
end

%% inverse planning for imrt
resultGUI  = matRad_fluenceOptimization(dij,cst,pln);

%% sequencing
if strcmp(pln.radiationMode,'photons') && (pln.propertiesOpt.runSequencing || pln.propertiesOpt.runDAO)
    %resultGUI = matRad_xiaLeafSequencing(resultGUI,stf,dij,5);
    %resultGUI = matRad_engelLeafSequencing(resultGUI,stf,dij,5);
    resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,5);
end

%% DAO
if strcmp(pln.radiationMode,'photons') && pln.propertiesOpt.runDAO
   resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);
   matRad_visApertureInfo(resultGUI.apertureInfo);
end

%% start gui for visualization of result
matRadGUI

%% indicator calculation and show DVH and QI
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI);

%% perform sampling
% select structures to include in sampling; leave empty to sample dose for all structures
structSel = {}; % structSel = {'PTV','OAR1'};
[caSampRes, mSampDose, plnSamp, resultGUInomScen] = matRad_sampling(ct,stf,cst,pln,resultGUI.w,structSel,[],[]);
[cstStat, resultGUIStat, param]                   = matRad_samplingAnalysis(ct,cst,plnSamp,caSampRes, mSampDose, resultGUInomScen,[]);



%% post processing
resultGUI = matRad_postprocessing(resultGUI, dij, pln, cst, stf);

%% export Plan
matRad_export_HITXMLPlan_modified('S02_orgpln',  pln, stf, resultGUI, 'stfMode')  %500000 minNbParticles HIT Minimum für Patienten, minNrParticlesIES, scan path mode: 'stfMode', 'backforth','TSP' (very slow)

%% calc 4D dose
[resultGUI, delivery] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI, 'S02_orgpln'); %'LiverDS221_1b_bf'); %'LiverDS221_1b_constRBE_bixel3_3_30fx60_bf'); %'Liver007_Bix5_1');  %'LiverDS221_1b_bf'); %'LiverDS221_wc5555_3mmBixel_bf'); %TKUH005_test');  
