% matRad example script 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_rc

% load patient data, i.e. ct, voi, cst
%load TG119.mat
load BOXPHANTOM.mat
%load LIVER.mat
%load PHANTOM_control.mat; ct.resolution.x = 2; ct.resolution.y = 2; ct.resolution.z = 2;

% meta information for treatment plan
pln.radiationMode   = 'protons';     % either photons / protons / carbon
%pln.machine         = 'generic_TOPAS_cropped';
pln.machine         = 'generic_MCsquare';


pln.numOfFractions  = 1;

% beam geometry settings
pln.propStf.bixelWidth              = 10; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longitudinalSpotSpacing = 10;
pln.propStf.gantryAngles            = 0; % [?] 
pln.propStf.couchAngles             = 0; % [?]
pln.propStf.numOfBeams              = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter               = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
%pln.propStf.isoCenter       = [51 0 51];
                            
% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]
%pln.propDoseCalc.doseGrid.resolution = ct.resolution;

%Turn on to correct for nozzle-to-skin air WEPL in analytical calculation
pln.propDoseCalc.airOffsetCorrection = true;

%Biology
modelName                   = 'none';
quantityOpt                 = 'physicalDose';  
pln.bioParam                = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% optimization settings
pln.propOpt.optimizer       = 'IPOPT';
                                      
                                                                           
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen'); % optimize on the nominal scenario     
 
pln.propMC.proton_engine = 'MCsquare';
%pln.propMC.proton_engine = 'TOPAS';


%% generate steering file
stf = matRad_generateStf(ct,cst,pln);
%stf = matRad_generateSingleBixelStf(ct,cst,pln); %Example to create a single beamlet stf

%% dose calculation
dij = matRad_calcParticleDose(ct, stf, pln, cst); %Calculate particle dose influence matrix (dij) with analytical algorithm
%dij = matRad_calcParticleDoseMC(ct,stf,pln,cst,1e4); %Calculate particle dose influence matrix (dij) with MC algorithm (slow!!)


resultGUI = matRad_fluenceOptimization(dij,cst,pln); %Optimize
%resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels,1),dij); %Use uniform weights

%resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
resultGUI_MC = matRad_calcDoseDirectMC(ct,stf,pln,cst,resultGUI.w,1e5);

%% Compare Dose
resultGUI = matRad_appendResultGUI(resultGUI,resultGUI_MC,0,'MC');
matRad_compareDose(resultGUI.physicalDose, resultGUI.physicalDose_MC, ct, cst, [1, 1, 0] , 'off', pln, [2, 2], 1, 'global');

%% GUI
matRadGUI
