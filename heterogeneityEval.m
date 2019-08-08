function [dose] = heterogeneityEval(particleType, doseMode, exportFolder, patients)
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

%%
% initialize matRad without clear command
addpath(genpath(pwd));

% set log level accordingly if you do _not_ want to do unit testing
if ~exist('unitTestBool','var') || ~unitTestBool
    
    param.calcDoseDirect = false;
    param.subIx          = [];
    param.logLevel       = 1;
    
    % set log level accordingly if want to do unit testing
else
    
    param.calcDoseDirect = false;
    param.subIx          = [];
    param.logLevel       = 3;
    
end

%%
p = genpath('../matRad');
addpath(p);

myFiles = dir(fullfile('dataImport/lowRes/*.mat'));
lungTissue={'Lung','GTV','PTV','CTV','ITV'};  % automatically assign HeterogeneityCorrection
% to the specified segmentations.
err = 0; % initialize Error counter
%for k = 1:length(myFiles)

prescribedDose = [60 30; 70 8; 66.41 6; 66.51 6]; % Prescribed doses and fractions for all 4 patients

for k = [patients]     % loop over all 4 patients
    
    %% Load in base data
    load(['dataImport/lowRes/',myFiles(k).name])
    
    pln.numOfFractions = prescribedDose(k,2);
    cst = matRad_slimCST(cst);
    cst = matRad_applyOARobjective(cst,prescribedDose(k));
    cstHetero = matRad_cstHeteroAutoassign(cst);
    
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
    
    if contains(doseMode,'Physical')
        modelName    = 'none';
        quantityOpt  = 'physicalDose';
        pln.propOpt.bioOptimization = 'none'; %'LEMIV_RBExD'; %'none'
    elseif contains(doseMode,'RBE')
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
    % try
    
    if contains(particleType,'HIT')        % HIT recalculated proton plan
        
        dose.hitHomo = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
        dose.hitHetero = matRad_calcDoseDirect(ct,stf,pln,cstHetero,resultGUI.w);
        
        mkdir(['dataImport/Results/',exportFolder,'/HIT/']);
        save(['dataImport/Results/',exportFolder,'/HIT/',myFiles(k).name(1:6),'.mat'])
        
    elseif contains(particleType,'Proton')        % reoptimized proton plan
        
        dijProt1 = matRad_calcParticleDose(ct,stf,pln,cst,param);
        dose.protHomo = matRad_fluenceOptimization(dijProt1,cst,pln,param);
        clear dijProt1
        dose.protHetero = matRad_calcDoseDirect(ct,stf,pln,cstHetero,dose.protHomo.w);
        
        mkdir(['dataImport/Results/',exportFolder,'/Protons/']);
        save(['dataImport/Results/',exportFolder,'/Protons/',myFiles(k).name(1:6),'.mat'])
        
    elseif contains(particleType,'Carbon')        % Carbon Ions
        
        pln.radiationMode = 'carbon';
        pln.machine = 'GenericAPM';
        pln.propStf.bixelWidth = 3;
        stfCarb = matRad_generateStf(ct,cst,pln);
        
        % Calculate dose for carbons
        dijCarb1 = matRad_calcParticleDose(ct,stfCarb,pln,cst);
        dose.carbHomo = matRad_fluenceOptimization(dijCarb1,cst,pln);
        clear dijCarb1
        dose.carbHetero = matRad_calcDoseDirect(ct,stfCarb,pln,cstHetero,dose.carbHomo.w);
        
        mkdir(['dataImport/Results/',exportFolder,'/Carbons/']);
        save(['dataImport/Results/',exportFolder,'/Carbons/',myFiles(k).name(1:6),'.mat'])
        
    else
        error('something is not right')
    end
    %  catch MException
    %       warning(['Error: Skipping Patient',myFiles(k).name(1:6),', Reason: ',MException.message]);
    %   end
end
%sendPushbullet('SUCCESS',['Simulation finished with ',num2str(err),' error(s)'])


end
