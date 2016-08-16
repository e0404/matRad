%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test COP approaches on different matrad patients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUT
clear all

numOfInputs    = 5;
%matRadPatients = {'TG119','PROSTATE','LIVER','BOXPHANTOM'};
matRadPatients = {'PROSTATE'};
root           = 'E:\Mescher\33_SourceShift_NumOfSparseSamplings';
resulSpec      = input('set result specification: ','s');

for m = 1:numOfInputs
    inputFolderTmp = strsplit(uigetdir(root),'\');
    inputFolder{m} = fullfile(inputFolderTmp{end-1},inputFolderTmp{end});

    % result folder
    resultDirTmp  = fullfile(root,'Results');
    foldername    = [datestr(now,'yyyy_mm_dd__HH_MM'),'_',inputFolderTmp{end},'_',resulSpec];
    resultDir{m}  = fullfile(resultDirTmp,foldername);

    if exist(resultDir{m},'dir') ~= 7
        mkdir(resultDir{m});
    end

    for i = 1:length(matRadPatients)
        if exist(fullfile(resultDir{m},matRadPatients{i}),'dir') ~= 7
        mkdir(fullfile(resultDir{m},matRadPatients{i}));
        end
    end
end

clearvars -except matRadPatients root resultDir inputFolder numOfInputs

%% Run Optimization
for m = 1:numOfInputs
    for i = 1:length(matRadPatients)
        % load matRad workspace
        matRadWorkspace = dir(fullfile(root,inputFolder{m},matRadPatients{i},'02_matRad_workspace','*.mat'));
        load(fullfile(root,inputFolder{m},matRadPatients{i},'02_matRad_workspace',matRadWorkspace.name))

        cstFiles = dir(fullfile(root,inputFolder{m},matRadPatients{i},'03_optimization_CST','*.mat'));
        for j = 1:length(cstFiles)
            % load cst
            load(fullfile(root,inputFolder{m},matRadPatients{i},'03_optimization_CST',cstFiles(j).name))

            % perform optimization
            resultGUI = matRad_fluenceOptimization(dij,cst,pln);
            save(fullfile(resultDir{m},matRadPatients{i},cstFiles(j).name),'ct','stf','pln','cst','resultGUI');
        end

    end
end

%% Run Validation
clear all 
addpath('E:\Mescher\33_SourceShift_NumOfSparseSamplings')

root      = 'E:\Mescher\33_SourceShift_NumOfSparseSamplings';
resultDir = dir(fullfile(root,'Results'));
resultDir = resultDir(cellfun(@(x) ~isequal(x,'.') & ~isequal(x,'..') ,{resultDir.name}));
resultDir = fullfile(root,'Results',{resultDir.name});

% load validation multScen structure
load(fullfile(root,'Input','ValidationMultScen_100Scenarios'))

for m = 1:length(resultDir)
    matRadPatients = dir(resultDir{m});
    matRadPatients = matRadPatients(cellfun(@(x) ~isequal(x,'.') & ~isequal(x,'..'),{matRadPatients.name}));
    matRadPatients = {matRadPatients.name};
    
    for i = 1:length(matRadPatients)
        optResults = dir(fullfile(resultDir{m},matRadPatients{i},'*.mat'));
        for j = 1:length(optResults)
            % load optimization result
            load(fullfile(resultDir{m},matRadPatients{i},optResults(j).name))
            
            % recalc dose with optimized w but with multScen structure
            [~,dijReCalc] = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w,multScen);
            
            % save recalulcated dose
            save(fullfile(resultDir{m},matRadPatients{i},[optResults(j).name(1:end-4),'_recalcDose.mat']),'dijReCalc');
                      
        end
    end
end

%% Plot validated and optimized DCH

clear all 
addpath('E:\Mescher\33_SourceShift_NumOfSparseSamplings')

root      = 'E:\Mescher\33_SourceShift_NumOfSparseSamplings';
resultDir = dir(fullfile(root,'Results'));
resultDir = resultDir(cellfun(@(x) ~isequal(x,'.') & ~isequal(x,'..') ,{resultDir.name}));
resultDir = fullfile(root,'Results',{resultDir.name});

% load validation multScen structure
load(fullfile(root,'Input','ValidationMultScen_100Scenarios'))

for m = 1:length(resultDir)
    matRadPatients = dir(resultDir{m});
    matRadPatients = matRadPatients(cellfun(@(x) ~isequal(x,'.') & ~isequal(x,'..'),{matRadPatients.name}));
    matRadPatients = {matRadPatients.name};
    
    for i = 1:length(matRadPatients)
        Results = dir(fullfile(resultDir{m},matRadPatients{i},'*.mat'));
        optResults = Results(cellfun(@(x) isempty(strfind(x,'recalcDose')),{Results.name}));
        valResults = Results(cellfun(@(x) ~isempty(strfind(x,'recalcDose')),{Results.name}));
        for j = 1:length(optResults)
            % plot validated DCH
            h = figure
            hold on
            load(fullfile(resultDir{m},matRadPatients{i},optResults(j).name))
            load(fullfile(resultDir{m},matRadPatients{i},valResults(j).name))
            for k = 1:length(dijReCalc.physicalDose)
                doseVec{k} = dijReCalc.physicalDose{k}(:,1);
            end
            plot_DCH({'prostate bed'},0.95,cst,doseVec,dijReCalc,'k');
            
            % plot optimized DCH
            Inputfolder     = strsplit(resultDir{m},'\');
            Inputfolder     = Inputfolder{end};
            Inputfolder     = strsplit(Inputfolder,'_');
            Inputfolder     = [Inputfolder{6},'_',Inputfolder{7}];
            matRadWorkspace = dir(fullfile(root,'Input',Inputfolder,matRadPatients{i},'02_matRad_workspace','*.mat'));
            
            load(fullfile(root,'Input',Inputfolder,matRadPatients{i},'02_matRad_workspace',matRadWorkspace.name))
            load(fullfile(resultDir{m},matRadPatients{i},optResults(j).name))
            
            doseVec = matRad_backProjection(resultGUI.w,dij,pln.bioOptimization);
            
            plot_DCH({'prostate bed'},0.95,cst,doseVec,dij,'b');
            
            legend({'validated DCH','DCH reference','optimized DCH'},'Location','southwest')
            
            savefig(h,fullfile(resultDir{m},matRadPatients{i},[optResults(j).name(1:end-4),'_comparison.fig']))
        end
    end
end
