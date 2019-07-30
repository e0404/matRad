% Heterogeneity correction evaluation script
%
% Computes the dose for the HIT plan (recalc), for protons and for carbons
% with and without heterogeneity correction. Saves the workspace for each
% calculation as .mat file.
%
% Set mode parameter to PhysicalDose or RBE for the calculation and correct
% folder assignment
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team.
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
p = genpath('../matRad');
addpath(p);
myFiles = dir(fullfile('dataImport/lowRes/*.mat'));
lungTissue={'Lung','GTV','PTV','CTV','ITV'};  % automatically assign HeterogeneityCorrection
% to the specified segmentations.
err = 0; % initialize Error counter
%for k = 1:length(myFiles)

mode = 'RBE';
dose = [60 30; 70 8; 66.41 6; 66.51 6]; % Prescribed doses and fractions for all 4 patients

for k = 1     % loop over all 4 patients
    
    %% Load in base data
    load(['dataImport/lowRes/',myFiles(k).name])
    
    pln.numOfFractions = dose(k,2);
    cst = matRad_slimCST(cst);
    cst = matRad_applyOARobjective(cst,dose(k));
    cstHetero = cst;
    
    lungIx =[];
    for i = 1:length(cst)
        if contains(cst{i,2},lungTissue)
            lungIx = [lungIx, i];
        end
    end
    
    for i = lungIx
        cstHetero{i,5}.HeterogeneityCorrection = 'Lung';
    end
    
    
    pln.propStf.bixelWidth = 3;
    
    %%
    % Define the biological optimization model for treatment planning along
    % with the quantity that should be used for optimization. Possible model values
    % are:
    % 'none':     physical optimization;
    % 'constRBE': constant RBE of 1.1;
    % 'MCN':      McNamara-variable RBE model for protons;
    % 'WED':      Wedenberg-variable RBE model for protons
    % 'LEM':      Local Effect Model
    % and possible quantityOpt are 'physicalDose', 'effect' or 'RBExD'.
    % As we use protons, we use a constant RBE of 1.1.
    
    if contains(mode,'Physical')
        modelName    = 'none';
        quantityOpt  = 'physicalDose';
        pln.propOpt.bioOptimization = 'none'; %'LEMIV_RBExD'; %'none'
    elseif contains(mode,'RBE')
        modelName    = 'LEM';
        quantityOpt  = 'RBExD';
        pln.propOpt.bioOptimization = 'LEMIV_RBExD'; %'LEMIV_RBExD'; %'none'
    else
        error(['Must choose valid dose cube']);
    end
    
    pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);
    pln.multScen = matRad_multScen(ct,'wcScen');
    
    pln.robOpt   = false;
    pln.sampling = false;
    
    %% Protons
    %try
%         %% HIT Plan protons
%         matRad_dispToConsole('Starting HIT calculation. \n',param,'info');
%         hitHomo = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
%         %%%
%         hitHetero = matRad_calcDoseDirect(ct,stf,pln,cstHetero,resultGUI.w);
%         
%         save(['dataImport/Results_',mode,'/HIT/',myFiles(k).name(1:6),'.mat'])
%         clearvars -except pln cst cstHetero ct lungTissue myFiles k stf param err
%         
%         %% Protons reoptimize
%         matRad_dispToConsole('Starting proton calculation. \n',param,'info');
%         dijProt1 = matRad_calcParticleDose(ct,stf,pln,cst,param);
%         protHomo = matRad_fluenceOptimization(dijProt1,cst,pln,param);
%         clear dijProt1
%         %%%
%         protHetero = matRad_calcDoseDirect(ct,stf,pln,cstHetero,protHomo.w);
%         
%         save(['dataImport/Results_',mode,'/Protons/',myFiles(k).name(1:6),'.mat'])
%         clearvars -except pln cst cstHetero ct lungTissue myFiles k stf param err
        
        %% Carbon Ions
        %%% Change to carbons
        
        pln.radiationMode = 'carbon';
        pln.machine = 'GenericAPM';
        pln.propStf.bixelWidth = 3;
        stfCarb = matRad_generateStf(ct,cst,pln);
        
        %% Calculate dose for carbons
        matRad_dispToConsole('Starting carbon calculation. \n',param,'info');
        dijCarb1 = matRad_calcParticleDose(ct,stfCarb,pln,cst);
        carbHomo = matRad_fluenceOptimization(dijCarb1,cst,pln);
        clear dijCarb1
        
        carbHetero = matRad_calcDoseDirect(ct,stfCarb,pln,cstHetero,carbHomo.w);
        
        save(['dataImport/Results_',mode,'/Carbons/',myFiles(k).name(1:6),'.mat'])
        sendPushbullet('SUCCESS',['Calculation complete for Patient',myFiles(k).name(1:6)])
        
%    catch MException
%        warning(['Error: Skipping Patient',myFiles(k).name(1:6),', Reason: ',MException.message]);
%        sendPushbullet('ERROR',['Skipping Patient',myFiles(k).name(1:6),', Reason: ',MException.message])
%       err = err + 1;
%    end
    close all
    clearvars -except lungTissue myFiles k param err dose
end
%sendPushbullet('SUCCESS',['Simulation finished with ',num2str(err),' error(s)'])



