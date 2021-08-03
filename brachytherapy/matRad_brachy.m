% matRad brachytherapy script
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

addpath(fileparts(pwd));
matRad_rc

% load patient data, i.e. ct, voi, cst
load PROSTATE.mat
matRad_changeCST

%% meta information for treatment plan
pln.numOfFractions                   = 1;

pln.radiationMode   = 'brachy';   
pln.machine         = 'HDR';        % 'LDR' or 'HDR' for brachy therapy
pln.propDoseCalc.TG43approximation = '2D'; %'1D' or '2D'
% make sure to always use HDR+2D or LDR+1D for the given base data


% needle settings
pln.propStf.needle.seedDistance      = 10; % [mm] seed distance on needle
pln.propStf.needle.seedsNo           = 6; % number of seeds per needle

%template grid definition (standard 13 x 13 grid with individual distances)
pln.propStf.template.activeNeedles = [0 0 0 1 0 1 0 1 0 1 0 0 0;... % 7.0
                                      0 0 1 0 1 0 0 0 1 0 1 0 0;... % 6.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 6.0
                                      1 0 1 0 1 0 0 0 1 0 1 0 1;... % 5.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 5.0
                                      1 0 1 0 1 0 0 0 1 0 1 0 1;... % 4.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 4.0
                                      1 0 1 0 1 0 0 0 1 0 1 0 1;... % 4.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 3.0
                                      1 0 1 0 1 0 1 0 1 0 1 0 1;... % 2.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 2.0
                                      1 0 1 0 1 0 0 0 0 0 1 0 1;... % 1.5
                                      0 0 0 0 0 0 0 0 0 0 0 0 0];   % 1.0
                                     %A a B b C c D d E e F f G

        

pln.propStf.template.normal      = [0,0,1];
pln.propStf.bixelWidth   = 5; % [mm] / also corresponds to
                              % template grid distance for brachy therapy
pln.propStf.isoCenter    = matRad_getIsoCenter(cst,ct,0); %  target center
pln.propStf.templateRoot = matRad_getTemplateRoot(ct,cst); % mass center of
%target in x and y and bottom in z

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

pln.propDoseCalc.DistanceCutoff    = 130; %[mm] % sets the maximum distance
% to which dose is calculated. 

% optimization settings
pln.propOpt.optimizer       = 'IPOPT';
pln.propOpt.bioOptimization = 'none';
        
                                              
% bookkeeping to ensure that GUI runs properly
pln.propOpt.runDAO          = false;  
pln.propOpt.runSequencing   = false; 
pln.propStf.gantryAngles    = [0:72:359]; 
pln.propStf.couchAngles     = [0 0 0 0 0]; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);

%% initial visualization and change objective function settings if desired
matRadGUI

%% generate steering file
stf = matRad_generateStf(ct,cst,pln,1);


%% dose calculation
dij = matRad_calcBrachyDose(ct,stf,pln,cst);


%% inverse planning for imrt or brachytherapy
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% start gui for visualization of result
matRadGUI

%% indicator calculation and show DVH and QI
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI);
