%% Example: Standrad Proton Treatment Plan recalculated with fine sampling 
%           and Monte Carlo
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
% (iii) how to recalculate the plan using the fine sampling algorithm
% (iv) how to recalculate the plan using MCsquare


%% Patient Data Import
% Let's begin with a clear Matlab environment and import the TG119
% patient into your workspace
clc,clear,close all;
load('TG119.mat');

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines the most 
% important cornerstones of your treatment plan.

%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon. In this
% example we would like to use protons for treatment planning. Next, we
% need to define a treatment machine to correctly load the corresponding 
% base data. matRad features generic base data in the file
% 'proton_Generic.mat'; consequently the machine has to be set accordingly
pln.radiationMode = 'protons';        
pln.machine       = 'generic_MCsquare';




%%
% The remaining plan parameters are set similarly to the previous example files
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [0, 45, 90];
pln.propStf.couchAngles   = [0,  0,  0];
pln.propStf.bixelWidth    = 50;
pln.propStf.longitudinalSpotSpacing = 100; 
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

pln.propOpt.optimizer       = 'IPOPT';
pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);


%% dose calculation
 % analytical dose without fine sampling
    dij = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI = matRad_calcCubes(ones(sum([stf(:).totalNumOfBixels]),1),dij);
    anaDose     = resultGUI.physicalDose;

%  % analytical dose with fine sampling
    % pln.propDoseCalc.fineSampling stores parameters defining the fine 
    % sampling simulation
       
    pln.propDoseCalc.fineSampling.method = 'russo'; 
    % method for weight calculation, availabe methods:
    %   'russo'
    %   'fitCircle', supports N = 2,3 and 8
    %   'fitSquare', supports N = 2 and 3
    
    pln.propDoseCalc.fineSampling.N = 21;
    % parameter to modify number of calculated FS sub beams
    %   'russo',        total number of beams = N^2
    %   'fitCircle',    total number of beams = (2*N + 1)^2
    %   'fitSquare',    total number of beams = (2^N - 1) * 6 + 1
    
    pln.propDoseCalc.fineSampling.sigmaSub = 2;
    % Gaussian standard deviation of the sub Gaussian beams, only used
    % when fine sampling method 'russo' is selected'
    
    dijFS = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI_FS = matRad_calcCubes(ones(sum([stf(:).totalNumOfBixels]),1),dijFS);
    resultGUI.physicalDoseFS = resultGUI_FS.physicalDose;
    anaFsDose   = resultGUI.physicalDoseFS;

 % Monte Carlo dose
    % Monte Carlo dose is set up to use base data set BDL_matRad.txt, which
    % is the data set, protons_genericMCsquare was fitted to
    resultGUI_MC = matRad_calcDoseDirectMC(ct,stf,pln,cst,ones(sum([stf(:).totalNumOfBixels]),1), 1e5);
    resultGUI.physicalDoseMC = resultGUI_MC.physicalDose;
    mcDose      = resultGUI.physicalDoseMC;

 % call matRad GUI
    matRadGUI

    
