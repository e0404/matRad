% Heterogeneity correction plotting script
%
% Generates evaluation plots based on the data calculated in
% heterogeneityEval.m 
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

mode = 'PhysicalDose';  

myFiles3 = dir(fullfile(['dataImport/Results_',mode,'/Carbons/*.mat']));

for n = 1:length(myFiles3)
    disp(['HIT: Calculating DVH for patient ',myFiles3(n).name(1:6)]);
    load(['dataImport/Results_',mode,'/HIT/',myFiles3(n).name(1:6)]);
    
    cstIndex = [3 11; 4 14; 4 18; 3 14];  % Position of Lung and PTV in the cst
    cst2 = cst(cstIndex(n,:),1:end);
    
    if contains(mode,'Physical')
        dvh1 = matRad_calcDVH(cst2,hitHomo.physicalDose);
        dvh2 = matRad_calcDVH(cst2,hitHetero.physicalDose);
    elseif contains(mode,'RBE')
        dvh1 = matRad_calcDVH(cst2,hitHomo.RBExD);
        dvh2 = matRad_calcDVH(cst2,hitHetero.RBExD);
    else
        error(['Must choose valid dose cube']);
    end

    dvhWindow = max([dvh1(1).doseGrid dvh2(1).doseGrid]);
    
    disp(['Plotting DVH for patient ',myFiles3(n).name(1:6)]);
    f = figure('Renderer', 'painters', 'Position', [10 100 1000 700]);
    set(f,'Color',[1 1 1]);
    matRad_showDVH(dvh1,cst2,pln);
    hold on
    matRad_showDVH(dvh2,cst2,pln,2);
    xlim([0 dvhWindow*1.1])
    title(['HIT Plan: DVH for pat. ',myFiles3(n).name(1:6)])
    saveas(f,['dataImport/Results_',mode,'/',myFiles3(n).name(1:6),'HIT.png'])
    
    clearvars -except n myFiles3 cstIndex cst2 mode
    %%
    disp(['Protons: Calculating DVH for patient ',myFiles3(n).name(1:6)]);
    load(['dataImport/Results_',mode,'/Protons/',myFiles3(n).name(1:6)]);
    
    if contains(mode,'Physical')
        dvh1 = matRad_calcDVH(cst2,protHomo.physicalDose);
        dvh2 = matRad_calcDVH(cst2,protHetero.physicalDose);
    elseif contains(mode,'RBE')
        dvh1 = matRad_calcDVH(cst2,protHomo.RBExD);
        dvh2 = matRad_calcDVH(cst2,protHetero.RBExD);
    else
        error(['Must choose valid dose cube']);
    end
    
    dvhWindow = max([dvh1(1).doseGrid dvh2(1).doseGrid]);
    
    disp(['Plotting DVH for patient ',myFiles3(n).name(1:6)]);
    f = figure('Renderer', 'painters', 'Position', [10 100 1000 700]);
    set(f,'Color',[1 1 1]);
    matRad_showDVH(dvh1,cst2,pln);
    hold on
    matRad_showDVH(dvh2,cst2,pln,2);
    xlim([0 dvhWindow*1.1])
    title(['MatRad Proton Plan: DVH for pat. ',myFiles3(n).name(1:6)])
    saveas(f,['dataImport/Results_',mode,'/',myFiles3(n).name(1:6),'MatRadProtons.png'])
    
    clearvars -except n myFiles3 cstIndex cst2 mode
    
    %%
    disp(['Carbons: Calculating DVH for patient ',myFiles3(n).name(1:6)]);
    load(['dataImport/Results_',mode,'/Carbons/',myFiles3(n).name(1:6)]);
    
    if contains(mode,'Physical')
        dvh1 = matRad_calcDVH(cst2,carbHomo.physicalDose);
        dvh2 = matRad_calcDVH(cst2,carbHetero.physicalDose);
    elseif contains(mode,'RBE')
        dvh1 = matRad_calcDVH(cst2,carbHomo.RBExD);
        dvh2 = matRad_calcDVH(cst2,carbHetero.RBExD);
    else
        error(['Must choose valid dose cube']);
    end
    
    dvhWindow = max([dvh1(1).doseGrid dvh2(1).doseGrid]);
    
    disp(['Plotting DVH for patient ',myFiles3(n).name(1:6)]);
    f = figure('Renderer', 'painters', 'Position', [10 100 1000 700]);
    set(f,'Color',[1 1 1]);
    matRad_showDVH(dvh1,cst2,pln);
    hold on
    matRad_showDVH(dvh2,cst2,pln,2);
    xlim([0 dvhWindow*1.1])
    title(['MatRad Carbon Plan: DVH for pat. ',myFiles3(n).name(1:6)])
    saveas(f,['dataImport/Results_',mode,'/',myFiles3(n).name(1:6),'MatRadCarbons.png'])
    
    clearvars -except n myFiles3 cstIndex cst2 mode
    
end