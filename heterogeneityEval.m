%% Computes Proton and Carbon dose as well as differences and quality indicators for every patient
matRad_rc
p = genpath('../matRad');
addpath(p);
myFiles = dir(fullfile('dataImport/lowRes/*.mat'));
lungTissue={'Lung','GTV','PTV','CTV','ITV'};
for k = 1:length(myFiles)
    %for k = 2
    
    %% Load in base data
    load(['dataImport/lowRes/',myFiles(k).name])
    
    pln.numOfFractions = 10;
    cst = matRad_slimCST(cst);
    cst = matRad_applyOARobjective(cst);
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
    
    pln.propOpt.bioOptimization = 'none'; %'LEMIV_RBExD'; %'none'
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
    
    modelName    = 'none';
    quantityOpt  = 'physicalDose';
    
    pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);
    pln.multScen = matRad_multScen(ct,'wcScen');
    
    pln.robOpt   = false;
    pln.sampling = false;
    
    %% Protons
  %  try
        %% HIT Plan protons
        matRad_dispToConsole('Starting HIT calculation. \n',param,'info');
      %  hitHomo = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
        %%%
        hitHetero = matRad_calcDoseDirect(ct,stf,pln,cstHetero,resultGUI.w);
        
        %protHomo = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
        %protHetero = matRad_calcDoseDirect(ct,stf,pln,cstHetero,resultGUI.w);
        
        %%% compare everything and plot
        %      close all
        %      [qI,gammaPassRate] = matRad_compareDose(hitHomo.physicalDose,hitHetero.physicalDose,ct,cst,[1 1 1 1],'off',pln);
        
        %%% Save all plots in subfolder
        %         FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
        %         for iFig = 1:length(FigList)
        %             FigHandle = FigList(iFig);
        %             FigName   = [myFiles(k).name(1:6),'_Fig',num2str(get(FigHandle, 'Number'))];
        %             savefig(FigHandle, ['dataImport/Protons/', FigName, '.fig']);
        %         end
        save(['dataImport/ResultWithPenalty/HIT/',myFiles(k).name(1:6),'.mat'])
        %      close all
        clearvars -except pln cst cstHetero ct lungTissue myFiles k stf param
        
        %% Protons reoptimize
        matRad_dispToConsole('Starting proton calculation. \n',param,'info');
        dijProt1 = matRad_calcParticleDose(ct,stf,pln,cst,param);
        protHomo = matRad_fluenceOptimization(dijProt1,cst,pln,param);
        clear dijProt1
        %%%
        protHetero = matRad_calcDoseDirect(ct,stf,pln,cstHetero,protHomo.w);
        
        %protHomo = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
        %protHetero = matRad_calcDoseDirect(ct,stf,pln,cstHetero,resultGUI.w);
        
        %%% compare everything and plot
        %       close all
        %      [qI,gammaPassRate] = matRad_compareDose(protHomo.physicalDose,protHetero.physicalDose,ct,cst,[1 1 1 1],'off',pln);
        
        %%% Save all plots in subfolder
        %         FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
        %         for iFig = 1:length(FigList)
        %             FigHandle = FigList(iFig);
        %             FigName   = [myFiles(k).name(1:6),'_Fig',num2str(get(FigHandle, 'Number'))];
        %             savefig(FigHandle, ['dataImport/Protons/', FigName, '.fig']);
        %         end
        save(['dataImport/ResultWithPenalty/Protons/',myFiles(k).name(1:6),'.mat'])
        %      close all
        clearvars -except pln cst cstHetero ct lungTissue myFiles k stf param
        
        %% Carbon Ions
        %%% Change to carbons
        
        pln.radiationMode = 'carbon';
        pln.machine = 'HIT_APM';
        pln.propOpt.bioOptimization = 'none'; %'LEMIV_RBExD'; %'none'
        pln.propStf.bixelWidth = 3;
        stfCarb = matRad_generateStf(ct,cst,pln);
        
        %% Calculate dose for carbons
        matRad_dispToConsole('Starting carbon calculation. \n',param,'info');
        dijCarb1 = matRad_calcParticleDose(ct,stfCarb,pln,cst);
        carbHomo = matRad_fluenceOptimization(dijCarb1,cst,pln);
        clear dijCarb1
        %%%
        %dijCarb2 = matRad_calcParticleDose(ct,stfCarb,pln,cstHetero);
        %carbHetero = matRad_fluenceOptimization(dijCarb2,cstHetero,pln);
        carbHetero = matRad_calcDoseDirect(ct,stfCarb,pln,cstHetero,carbHomo.w);
        
        %%% recalc plan 1
        %carbHetereo = matRad_calcCubes(carbHomo.w,dijCarb2,cstHetero);
        
        %%%
        %      close all
        %      [qI,gammaPassRate]=matRad_compareDose(carbHomo.physicalDose,carbHetero.physicalDose,ct,cst,[1 1 1 1],'off',pln);
        %      matRad_plotSOBP(ct,stfCarb,carbHomo.physicalDose,carbHetero.physicalDose)
        
        %%% Save all plots in subfolder
        %         FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
        %         for iFig = 1:length(FigList)
        %             FigHandle = FigList(iFig);
        %             FigName   = [myFiles(k).name(1:6),'_Fig',num2str(get(FigHandle, 'Number'))];
        %             savefig(FigHandle, ['dataImport/Carbons/', FigName, '.fig']);
        %         end
        save(['dataImport/ResultWithPenalty/Carbons/',myFiles(k).name(1:6),'.mat'])
        sendPushbullet('SUCCESS',['Calculation complete for Patient',myFiles(k).name(1:6)])
        
 %   catch MException
%        warning(['Error: Skipping Patient',myFiles(k).name(1:6),', Reason: ',MException.message]);
 %       sendPushbullet('ERROR',['Skipping Patient',myFiles(k).name(1:6),', Reason: ',MException.message])
%    end
    %  close all
    clearvars -except lungTissue myFiles k param
    
end
