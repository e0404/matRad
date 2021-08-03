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

%pln.radiationMode   = 'brachy';     % either 'photons' / 'protons' / 'carbon' / 'brachy'
%pln.machine         = 'HDR';        % 'Generic' for photons, protons, carbon // 'LDR' or 'HDR' for brachy
pln.radiationMode   = 'photons';     % either photons / protons / carbon
pln.machine         = 'Generic';

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

% optimization settings
pln.propOpt.optimizer       = 'IPOPT';
pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                              % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose

% beam specific settings (ignore for brachy)
        pln.numOfFractions  = 30;

        % beam geometry settings
        pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
        pln.propStf.gantryAngles    = [0:72:359]; % [?]
        pln.propStf.couchAngles     = [0 0 0 0 0]; % [?]
        pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
        pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
        % photon optimization
        pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
        pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

% brachy specific settings (ignore for external beam)        
        pln.numOfFractions                   = 1;
        
        % geometry settings
        pln.propStf.templateRoot             = matRad_getTemplateRoot(ct,cst); % mass center of target in x and y and bottom in z
        pln.propStf.needle.seedDistance      = 10; % [mm] seed distance on needle
        pln.propStf.needle.seedsNo           = 11; % number of seeds per needle
        % Option one: template grid definition
            pln.propStf.template.numOfXPoints  = 8;
            pln.propStf.template.numOfYPoints = 6;
            pln.propStf.template.xScale        = 25; % [mm] distance of neighbouring points
            pln.propStf.template.yScale       = 25; % [mm] distance of neighbouring points
        % Option two: direct array of template points
        % put in eather 'none' or 4xN array, where each column is the
        % 3+1 vector of one template hole (fourth entry always one, z coord - third entry zero)
            pln.propStf.templateDirect          = 'none';

        % displacement - rotation matrix
            % [4x4] matrix. First three columns:
            % orthonormal base vectors of rotated template and zero in the fourth entry
            % Fourth column: template origin in LPS mm coordinate system
            % and one in the fourth entry
        pln.propStf.shiftRotMtx = [1 0 0 pln.propStf.templateRoot(1);...
                                  0 1 0 pln.propStf.templateRoot(2);...
                                  0 0 1 pln.propStf.templateRoot(3);...
                                  0 0 0 1];
        
        % dose calculation settings
        pln.propDoseCalc.durationImplanted = Inf; % only important LDR therapy
        pln.propDoseCalc.TG43approximation = '2D'; %'1D' or '2D' 
        
%% initial visualization and change objective function settings if desired
matRadGUI

%% generate steering file
switch pln.radiationMode
    case {'photons','protons','carbon'}
stf = matRad_generateStf(ct,cst,pln);
    case 'brachy'
stf = matRad_generateStf(ct,cst,pln);
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

%% inverse planning for imrt or brachytherapy
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

