myFiles = dir(fullfile('dataImport/ResultWithPenalty/HIT/*.mat'));


for n = 1:length(myFiles)
    disp(['HIT: Calculating DVH for patient ',myFiles(n).name(1:6)]);
    load(['dataImport/ResultWithPenalty/HIT/',myFiles(n).name(1:6)]);
    
    cstIndex = [3 11; 4 14; 4 18; 4 13; 3 14];
    cst2 = cst(cstIndex(n,:),1:end);
    
    dvh1 = matRad_calcDVH(cst2,hitHomo.physicalDose/2);
    dvh2 = matRad_calcDVH(cst2,hitHetero.physicalDose/2);
    pln.bioParam.model = 'none';
    
    dvhWindow = max([dvh1(1).doseGrid dvh2(1).doseGrid]);
    
    disp(['Plotting DVH for patient ',myFiles(n).name(1:6)]);
    f = figure('Renderer', 'painters', 'Position', [10 100 1000 700]);
    set(f,'Color',[1 1 1]);
    matRad_showDVH(dvh1,cst2,pln);
    hold on
    matRad_showDVH(dvh2,cst2,pln,2);
    xlim([0 dvhWindow*1.1])
    title(['HIT Plan: DVH for pat. ',myFiles(n).name(1:6)])
    saveas(f,['dataImport/ResultWithPenalty/',myFiles(n).name(1:6),'HIT.png'])
    
    clearvars -except n myFiles cstIndex cst2
%%
    disp(['Protons: Calculating DVH for patient ',myFiles(n).name(1:6)]);
    load(['dataImport/ResultWithPenalty/Protons/',myFiles(n).name(1:6)]);
      
    dvh1 = matRad_calcDVH(cst2,protHomo.physicalDose);
    dvh2 = matRad_calcDVH(cst2,protHetero.physicalDose);
    pln.bioParam.model = 'none';
    
    dvhWindow = max([dvh1(1).doseGrid dvh2(1).doseGrid]);
    
    disp(['Plotting DVH for patient ',myFiles(n).name(1:6)]);
    f = figure('Renderer', 'painters', 'Position', [10 100 1000 700]);
    set(f,'Color',[1 1 1]);
    matRad_showDVH(dvh1,cst2,pln);
    hold on
    matRad_showDVH(dvh2,cst2,pln,2);
    xlim([0 dvhWindow*1.1])
    title(['MatRad Proton Plan: DVH for pat. ',myFiles(n).name(1:6)])
    saveas(f,['dataImport/ResultWithPenalty/',myFiles(n).name(1:6),'MatRadProtons.png'])
    
    clearvars -except n myFiles cstIndex cst2

%%
    disp(['Carbons: Calculating DVH for patient ',myFiles(n).name(1:6)]);
    load(['dataImport/ResultWithPenalty/Carbons/',myFiles(n).name(1:6)]);
    
    dvh1 = matRad_calcDVH(cst2,carbHomo.physicalDose);
    dvh2 = matRad_calcDVH(cst2,carbHetero.physicalDose);
    pln.bioParam.model = 'none';
    
    dvhWindow = max([dvh1(1).doseGrid dvh2(1).doseGrid]);
    
    disp(['Plotting DVH for patient ',myFiles(n).name(1:6)]);
    f = figure('Renderer', 'painters', 'Position', [10 100 1000 700]);
    set(f,'Color',[1 1 1]);
    matRad_showDVH(dvh1,cst2,pln);
    hold on
    matRad_showDVH(dvh2,cst2,pln,2);
    xlim([0 dvhWindow*1.1])
    title(['MatRad Carbon Plan: DVH for pat. ',myFiles(n).name(1:6)])
    saveas(f,['dataImport/ResultWithPenalty/',myFiles(n).name(1:6),'MatRadCarbons.png'])

    clearvars -except n myFiles cstIndex cst2

end