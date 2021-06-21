% matRad script
%
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
matRad_rc

% load patient data, i.e. ct, voi, cst

%load HEAD_AND_NECK
%load TG119.mat
load PROSTATE.mat
%load LIVER.mat
%load BOXPHANTOM.mat

%% meta information for treatment plan

pln.radiationMode   = 'brachy';     % either photons / protons / carbon / brachy
%pln.radiationMode   = 'protons';  
pln.machine         = 'LDR';

switch pln.radiationMode
    case {'photons','protons','carbon'}
        
        pln.numOfFractions  = 30;
        
        % beam geometry settings
        pln.propStf.bixelWidth       = 5; % [mm] / also corresponds to lateral spot spacing for particles
        pln.propStf.gantryAngles     = [0:72:359]; % [?]
        pln.propStf.couchAngles      = [0 0 0 0 0]; % [?]
        pln.propStf.numOfBeams       = numel(pln.propStf.gantryAngles);
        pln.propStf.isoCenter        = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
        
        % dose calculation settings
        pln.propDoseCalc.doseGrid.resolution.x = 10; % [mm]
        pln.propDoseCalc.doseGrid.resolution.y = 10; % [mm]
        pln.propDoseCalc.doseGrid.resolution.z = 10; % [mm]
        
        % optimization settings
        pln.propOpt.optimizer       = 'IPOPT';
        pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
        pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
        pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

    case 'brachy'
        % template geometry settings
        pln.numOfFractions                   = 1;
        pln.propStf.template.numOfHorPoints  = 8;
        pln.propStf.template.numOfVertPoints = 6;
        pln.propStf.template.Xscale          = 20; % [mm]
        pln.propStf.template.Yscale          = 20; % [mm]
        pln.propStf.needle.seedDistance      = 12; % [mm]
        pln.propStf.needle.seedsNo           = 11; 
            %unit vectors of displaced, rotated template coordinate system
        pln.propStf.orientation.Xdir = normalize([1,0,0],'norm');
        pln.propStf.orientation.Ydir = normalize([0,1,0],'norm');
        pln.propStf.orientation.Zdir = cross(pln.propStf.orientation.Xdir,pln.propStf.orientation.Ydir);
        pln.propStf.orientation.offset = [-83,-408,50]; % [mm]
        assert(pln.propStf.orientation.Xdir*pln.propStf.orientation.Ydir' == 0,'Xdir and Ydir are not orthogonal') %throw error if directions are not orthogonal
        
        % dose calculation settings
        pln.propDoseCalc.durationImplanted = Inf;
        pln.propDoseCalc.TG43approximation = '1D'; %'1D' or '2D'
        pln.propDoseCalc.doseGrid.resolution.x = 10; % [mm]
        pln.propDoseCalc.doseGrid.resolution.y = 10; % [mm]
        pln.propDoseCalc.doseGrid.resolution.z = 10; % [mm]
        
        % optimization settings
        pln.propOpt.optimizer           = 'IPOPT';
        pln.propOpt.bioOptimization     = 'none';
        pln.propOpt.runDAO              = false;
        pln.propOpt.runSequencing       = false;
        pln.propOpt.lowerWeightBounds   = 0;
        pln.propOpt.upperWeightBounds   = 1;
end

%% initial visualization and change objective function settings if desired
matRadGUI

%% generate steering file
switch pln.radiationMode
    case {'photons','protons','carbon'}
        stf = matRad_generateStf(ct,cst,pln);
    case 'brachy'
        stf = matRad_generateBrachyStf(ct,cst,pln);
end


%% dose calculation
switch pln.radiationMode
    case 'photons'
    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
    %dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
    case {'protons','carbon'}
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
    case 'brachy'
    dij = matRad_calcBrachyDose(ct,stf,pln,cst);
end

%% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% sequencing
if strcmp(pln.radiationMode,'photons') && (pln.propOpt.runSequencing || pln.propOpt.runDAO)
    %resultGUI = matRad_xiaLeafSequencing(resultGUI,stf,dij,5);
    %resultGUI = matRad_engelLeafSequencing(resultGUI,stf,dij,5);
    resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,5);
end

%% DAO
if strcmp(pln.radiationMode,'photons') && pln.propOpt.runDAO
   resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);
   matRad_visApertureInfo(resultGUI.apertureInfo);
end

%% start gui for visualization of result
matRadGUI

%% indicator calculation and show DVH and QI
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI);

