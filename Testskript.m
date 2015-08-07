% % TestScript

%% execute matRad until Sequencing
    clear
     close all
     clc

    % load patient data, i.e. ct, voi, cst

    load HEAD_AND_NECK
    %load TG119.mat
    %load PROSTATE.mat
%     load LIVER.mat
    %load BOXPHANTOM.mat

    % meta information for treatment plan
    pln.SAD             = 1000; %[mm]
    pln.isoCenter       = matRad_getIsoCenter(cst,ct,0);
    pln.bixelWidth      = 15; % [mm] / also corresponds to lateral spot spacing for particles
    pln.gantryAngles    = [0:90:359]; % [°]
    pln.couchAngles     = [0:90:359]; % [°]
    pln.numOfBeams      = numel(pln.gantryAngles);
    pln.numOfVoxels     = numel(ct.cube);
    pln.voxelDimensions = size(ct.cube);
    pln.radiationMode   = 'photons'; % either photons / protons / carbon
    pln.bioOptimization = 'none'; % none: physical optimization; effect: effect-based optimization; RBExD: optimization of RBE-weighted dose
    pln.numOfFractions  = 1;

    %% initial visualization and change objective function settings if desired
% skip the GUI
%     matRadGUI

    %% generate steering file
    stf = matRad_generateStf(ct,cst,pln);

    %% dose calculation
    if strcmp(pln.radiationMode,'photons')
        dij = matRad_calcPhotonDose(ct,stf,pln,cst,0);
    elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
        dij = matRad_calcParticleDose(ct,stf,pln,cst,0);
    end

    %% inverse planning for imrt
    resultGUI = matRad_fluenceOptimization(dij,cst,pln,1);

    %% sequencing
    if strcmp(pln.radiationMode,'photons')
%         Sequencing = matRad_xiaLeafSequencing(resultGUI.w,stf,7,1);
        Sequencing = matRad_engelLeafSequencing(resultGUI.w,stf,7);
        resultGUI = matRad_mxCalcDose(dij,Sequencing.w,cst);
    end
    
%% get information from sequencing and visualize

shapeInfo = tk_getParameters(Sequencing,stf,pln,0);
tk_visualizeMLC(shapeInfo,pln)

% shapeInfo = tk_getSequencingParameters(Sequencing,pln,stf,1);
% tk_drawShapes(shapeInfo,pln);

