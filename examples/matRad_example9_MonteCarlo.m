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

% meta information for treatment plan
pln.radiationMode   = 'protons';     % either photons / protons / carbon

pln.machine          = 'Generic';
%pln.machine          = 'generic_MCsquare';


pln.numOfFractions  = 1;

% beam geometry settings
pln.propStf.bixelWidth              = 10; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longitudinalSpotSpacing = 10;
pln.propStf.gantryAngles            = 0; % [?] 
pln.propStf.couchAngles             = 0; % [?]
pln.propStf.numOfBeams              = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter               = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
%pln.propStf.isoCenter              = [51 0 51];
                            
% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]
%pln.propDoseCalc.doseGrid.resolution = ct.resolution;

%Activate/Deactive the minimal correction for air WEPL given the SSD 
%differences between base data fit and treatment plan
%
%pln.propDoseCalc.airOffsetCorrection = false;

% optimization settings
pln.propOpt.optimizer       = 'IPOPT';
pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

%ln.propMC.proton_engine    = 'TOPAS'; %Requires separate topas installation
pln.propMC.proton_engine    = 'MCsquare';

%Enable LET calculation - Does not work yet for MCsquare data
pln.propDoseCalc.calcLET    = false;

%Enable/Disable use of range shifter (has effect only when we need to fill 
%up the low-range region)
pln.propStf.useRangeShifter = true;  

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);
%stf = matRad_generateSingleBixelStf(ct,cst,pln); %Example to create a single beamlet stf

%% dose calculation

dij = matRad_calcParticleDose(ct, stf, pln, cst); %Calculate particle dose influence matrix (dij) with analytical algorithm
%dij = matRad_calcParticleDoseMC(ct,stf,pln,cst,1e4); %Calculate particle dose influence matrix (dij) with MC algorithm (slow!!)


resultGUI = matRad_fluenceOptimization(dij,cst,pln); %Optimize
%resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels,1),dij); %Use uniform weights


%% MC calculation
%resultGUI_recalc = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);       %Recalculate particle dose analytically
resultGUI_recalc = matRad_calcDoseDirectMC(ct,stf,pln,cst,resultGUI.w,1e6); %Recalculate particle dose with MC algorithm
resultGUI.physicalDose_MC = resultGUI_recalc.physicalDose;

%% Compare Dose
matRad_compareDose(resultGUI.physicalDose, resultGUI_recalc.physicalDose, ct, cst, [1, 1, 0] , 'off', pln, [2, 2], 1, 'global');

%% Compare LET
if isfield(resultGUI,'LET') && isfield(resultGUI_recalc,'LET')
    matRad_compareDose(resultGUI.LET, resultGUI_recalc.LET, ct, cst, [1, 1, 0] , 'off', pln, [2, 2], 1, 'global');
    resultGUI.LET_MC = resultGUI_recalc.LET;
end

%% GUI
matRadGUI


