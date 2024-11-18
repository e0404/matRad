% matRad example script 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
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
pln.machine         = 'Generic';

%Biology
pln.bioModel = 'none';
%Scenario Model
pln.multScen = 'nomScen';

pln.numOfFractions  = 1;

% beam geometry settings
pln.propStf.bixelWidth              = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longitudinalSpotSpacing = 3;
pln.propStf.gantryAngles            = 0; % [?] 
pln.propStf.couchAngles             = 0; % [?]
pln.propStf.numOfBeams              = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter               = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
%pln.propStf.isoCenter       = [51 0 51];
                            
% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]
%pln.propDoseCalc.doseGrid.resolution = ct.resolution;

%Turn on to correct for nozzle-to-skin air WEPL in analytical calculation
pln.propDoseCalc.airOffsetCorrection = true;

% optimization settings
pln.propOpt.optimizer       = 'IPOPT';
                                                                           
%Enable/Disable use of range shifter (has effect only when we need to fill 
%up the low-range region)
pln.propStf.useRangeShifter = false;  
pln.propStf.generator = 'ParticleSingleSpot'; 

%Enable LET calculation
pln.propDoseCalc.calcLET = true;
pln.propOpt.quantityOpt = 'physicalDose';

% Enable/Disable local computation with TOPAS. Enabling this will generate
% the necessary TOPAS files to run the simulation on any machine or server.
% pln.propMC.externalCalculation = true;

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);

%% analytical dose calculation
pln.propDoseCalc.engine = 'HongPB';
%pln.propDoseCalc.numHistoriesPerBeamlet = 1e6;

dij = matRad_calcDoseInfluence(ct, cst,stf, pln); %Calculate particle dose influence matrix (dij) with analytical algorithm

resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels,1),dij); %Use uniform weights
%resultGUI = matRad_fluenceOptimization(dij,cst,pln); %Optimize


%% Monte Carlo dose calculation
% select Monte Carlo engine ('MCsquare' very fast for physical protons, 'TOPAS' slow but versatile for everything else)
pln.propDoseCalc.engine = 'MCsquare';
%pln.propDoseCalc.engine = 'TOPAS';

% set number of histories lower than default for this example (default: 1e8)
pln.propDoseCalc.numHistoriesDirect = 5e6;
%pln.propDoseCalc.externalCalculation = 'write';
resultGUI_MC = matRad_calcDoseForward(ct,cst,stf,pln,resultGUI.w);

%% Read an external calculation with TOPAS if externalCalculation was set to 'write'
%pln.propDoseCalc.externalCalculation = resultGUI_MC.meta.TOPASworkingDir;
%resultGUI_MC = matRad_calcDoseForward(ct,cst,stf,pln,resultGUI.w);

%% Compare Dose
resultGUI = matRad_appendResultGUI(resultGUI,resultGUI_MC,true,pln.propDoseCalc.engine);
matRad_compareDose(resultGUI.physicalDose, resultGUI.(['physicalDose_' pln.propDoseCalc.engine]), ct, cst, [1, 1, 0] , 'off', pln, [2, 2], 3, 'global');


%% Compare LET
if isfield(resultGUI,'LET') && isfield(resultGUI_MC,'LET')
    matRad_compareDose(resultGUI.LET, resultGUI.(['LET_' pln.propDoseCalc.engine]), ct, cst, [1, 1, 0] , 'off', pln, [2, 2], 1, 'global');    
end

%% GUI
matRadGUI
