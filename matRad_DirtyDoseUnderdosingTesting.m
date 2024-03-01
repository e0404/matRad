%% Proton Optimization 2 beam without dirty Dose
matRad_rc;
clear("dij","pln","resultGUI","ct","cst","stf")
load("PROSTATE.mat")

cube = zeros(183,183,90);
cube(cst{6,4}{1}) = 1;
vResolution = ct.resolution;
vMargin = [];
vMargin.x = 5;
vMargin.y = 5;
vMargin.z = 5;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);

cst{11,1}    = 10;
cst{11,2}    = 'Margin';
cst{11,3}    = 'OAR';
cst{11,4}{1} = find(mVOIEnlarged);
cst{11,5}    = cst{3,5};

cst{11,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,60)); 

cube = zeros(183,183,90);
cube(cst{1,4}{1}) = 1;
vResolution = ct.resolution;
vMargin = [];
vMargin.x = 5;
vMargin.y = 5;
vMargin.z = 5;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);

cst{12,1}    = 11;
cst{12,2}    = 'Margin';
cst{12,3}    = 'OAR';
cst{12,4}{1} = find(mVOIEnlarged);
cst{12,5}    = cst{3,5};

cst{12,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,60)); 
%%
cst{8,6}{2} = struct(LETtObjectives.matRad_SquaredOverdosingLETt(100,4));
cst{6,6}{2} = struct(LETdObjectives.matRad_SquaredUnderdosingLETd(100,4));
cst{1,6}{2} = struct(LETdObjectives.matRad_SquaredOverdosingLETd(100,20));
% cst{6,6}{3} = struct(DirtyDoseObjectives.matRad_ClusterDirtyDoseVariance(100));
%%
%cst{6,6}{2} = struct(mLETDoseObjectives.matRad_SquaredUnderdosingmLETDose(100,10));


%% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;


% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

%% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);
dij = matRad_calcLETvD(dij);

% Inverse Optimization for IMPT
%%
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
dvh = matRad_calcDVH(cst,resultGUI.physicalDose);
matRad_indicatorWrapper(cst,pln,resultGUI)
%matRad_showDVH(dvh,cst);
%% DVH
fig = figure;
subplot(2,1,1)
for i = [1 6 8]
    plot(dvh_normalMphysDose(i).doseGrid,dvh_normalMphysDose(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','-','LineWidth',1.5)
    hold on
    plot(dvh_mLETDoseOverdosing_physDose(i).doseGrid,dvh_mLETDoseOverdosing_physDose(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','--','LineWidth',1.5)
    hold on
    plot(dvh_MDirtyDoseOverdosing_physDose(i).doseGrid,dvh_MDirtyDoseOverdosing_physDose(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle',':','LineWidth',1.5)
    % hold on
    % plot(dvh_mLETDoseOverdosing_mLETDose120_physDose(i).doseGrid,dvh_mLETDoseOverdosing_mLETDose120_physDose(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','-.','LineWidth',1.5)
    % hold on
    % plot(dvh_MDirtyDoseOverdosing_DirtyDose30_physDose(i).doseGrid,dvh_MDirtyDoseOverdosing_DirtyDose30_physDose(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','none','Marker',"x",'MarkerSize',3)
    % 
end
title('DVH')
xlabel('physical dose in Gy')
ylabel('Volume in %')
% legend('Rectum normal','Rectum LETd','Rectum dirty dose','Rectum LETd & Target LETd','Rectum dirtyDose & Target dirtyDose','PTV 68 normal','PTV 68 LETd', ...
%     'PTV 68 dirty dose','PTV 68 LETd & Target LETd','PTV 68 dirtyDose & Target dirtyDose','Bladder normal','Bladder LETd', ...
%     'Bladder dirty dose','Bladder LETd & Target LETd','Bladder dirtyDose & Target dirtyDose','Body normal','Body LETd','Body dirty dose', ...
%     'Body LETd & Target LETd','Body dirtyDose & Target dirtyDose',Location='eastoutside')

subplot(2,1,2)
for i = [1 6 8]
    plot(dvh_normalMdirtyDose(i).doseGrid,dvh_normalMdirtyDose(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','-','LineWidth',1.5)
    hold on
    plot(dvh_mLETDoseOverdosing_dirtyDose(i).doseGrid,dvh_mLETDoseOverdosing_dirtyDose(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','--','LineWidth',1.5)
    hold on
    plot(dvh_MDirtyDoseOverdosing_dirtyDose(i).doseGrid,dvh_MDirtyDoseOverdosing_dirtyDose(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle',':','LineWidth',1.5)
    % hold on
    % plot(dvh_mLETDoseOverdosing_mLETDose120_dirtyDose(i).doseGrid,dvh_mLETDoseOverdosing_mLETDose120_dirtyDose(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','-.','LineWidth',1.5)
    % hold on
    % plot(dvh_MDirtyDoseOverdosing_DirtyDose30_dirtyDose(i).doseGrid,dvh_MDirtyDoseOverdosing_DirtyDose30_dirtyDose(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','none','Marker',"x",'MarkerSize',3)

end

xlabel('DD in Gy')
xlim = 3;
ylabel('Volume in %')
% legend('Rectum normal','Rectum LETd','Rectum dirty dose','Rectum LETd & Target LETd','Rectum dirtyDose & Target dirtyDose','PTV 68 normal','PTV 68 LETd', ...
%     'PTV 68 dirty dose','PTV 68 LETd & Target LETd','PTV 68 dirtyDose & Target dirtyDose','Bladder normal','Bladder LETd', ...
%     'Bladder dirty dose','Bladder LETd & Target LETd','Bladder dirtyDose & Target dirtyDose','Body normal','Body LETd','Body dirty dose', ...
%     'Body LETd & Target LETd','Body dirtyDose & Target dirtyDose',Location='eastoutside')

% subplot(3,1,3)
% for i = [1 6 8 9]
%     plot(dvh_normalMLETd(i).doseGrid,dvh_normalMLETd(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','-','LineWidth',1.5)
%     hold on
%     plot(dvh_mLETDoseOverdosing_LETd(i).doseGrid,dvh_mLETDoseOverdosing_LETd(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','--','LineWidth',1.5)
%     hold on
%     plot(dvh_MDirtyDoseOverdosing_LETd(i).doseGrid,dvh_MDirtyDoseOverdosing_LETd(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle',':','LineWidth',1.5)
%     hold on
%     plot(dvh_mLETDoseOverdosing_mLETDose120_LETd(i).doseGrid,dvh_mLETDoseOverdosing_mLETDose120_LETd(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','-.','LineWidth',1.5)
%     hold on
%     plot(dvh_MDirtyDoseOverdosing_DirtyDose30_LETd(i).doseGrid,dvh_MDirtyDoseOverdosing_DirtyDose30_LETd(i).volumePoints,'Color',[cst{i,5}.visibleColor(1,1),cst{i,5}.visibleColor(1,2),cst{i,5}.visibleColor(1,3)],'LineStyle','none','Marker',"x",'MarkerSize',3)
% 
% end
% 
% xlim = 10;
% xlabel('LETd in keV/µm')
% ylabel('Volume in %')
% 
% fig = gcf;

% add legend
Lgnd = legend('Rectum Ref','Rectum LETxD','Rectum DD','PTV Ref','PTV LETxD', ...
    'PTV 68 DD','Bladder Ref','Bladder LETxD', ...
    'Bladder DD',Location='southoutside',NumColumns=1);

Lgnd.Position(1) = -0.012;
Lgnd.Position(2) = 0.01;

%% DVH
fig = figure;
linestyle = {'-','--','-.',':'};
subplot(2,1,1)
for i = [1 6 8]
    for j = [1 2 3 4]
        plot(dvh_normalMphysDose(i).doseGrid,dvh_normalMphysDose(i).volumePoints,'Color','red','LineStyle',linestyle(j),'LineWidth',1.5)
        hold on
        plot(dvh_mLETDoseOverdosing_physDose(i).doseGrid,dvh_mLETDoseOverdosing_physDose(i).volumePoints,'Color','blue','LineStyle',linestyle(j),'LineWidth',1.5)
        hold on
        plot(dvh_MDirtyDoseOverdosing_physDose(i).doseGrid,dvh_MDirtyDoseOverdosing_physDose(i).volumePoints,'Color','green','LineStyle',linestyle(j),'LineWidth',1.5)
        % hold on
        % plot(dvh_mLETDoseOverdosing_mLETDose120_physDose(i).doseGrid,dvh_mLETDoseOverdosing_mLETDose120_physDose(i).volumePoints,'Color','magenta','LineStyle',linestyle(j),'LineWidth',1.5)
        % hold on
        % plot(dvh_MDirtyDoseOverdosing_DirtyDose30_physDose(i).doseGrid,dvh_MDirtyDoseOverdosing_DirtyDose30_physDose(i).volumePoints,'Color','black','LineStyle',linestyle(j),'LineWidth',1.5)
    end
end
title('DVH')
xlabel('physical dose in Gy')
ylabel('Volume in %')
% legend('Rectum normal','Rectum LETd','Rectum dirty dose','Rectum LETd & Target LETd','Rectum dirtyDose & Target dirtyDose','PTV 68 normal','PTV 68 LETd', ...
%     'PTV 68 dirty dose','PTV 68 LETd & Target LETd','PTV 68 dirtyDose & Target dirtyDose','Bladder normal','Bladder LETd', ...
%     'Bladder dirty dose','Bladder LETd & Target LETd','Bladder dirtyDose & Target dirtyDose','Body normal','Body LETd','Body dirty dose', ...
%     'Body LETd & Target LETd','Body dirtyDose & Target dirtyDose',Location='eastoutside')

subplot(2,1,2)
for i = [1 6 8]
    for j = [1 2 3 4]
        plot(dvh_normalMdirtyDose(i).doseGrid,dvh_normalMdirtyDose(i).volumePoints,'Color','red','LineStyle',linestyle(j),'LineWidth',1.5)
        hold on
        plot(dvh_mLETDoseOverdosing_dirtyDose(i).doseGrid,dvh_mLETDoseOverdosing_dirtyDose(i).volumePoints,'Color','blue','LineStyle',linestyle(j),'LineWidth',1.5)
        hold on
        plot(dvh_MDirtyDoseOverdosing_dirtyDose(i).doseGrid,dvh_MDirtyDoseOverdosing_dirtyDose(i).volumePoints,'Color','green','LineStyle',linestyle(j),'LineWidth',1.5)
        % hold on
        % plot(dvh_mLETDoseOverdosing_mLETDose120_dirtyDose(i).doseGrid,dvh_mLETDoseOverdosing_mLETDose120_dirtyDose(i).volumePoints,'Color','magenta','LineStyle',linestyle(j),'LineWidth',1.5)
        % hold on
        % plot(dvh_MDirtyDoseOverdosing_DirtyDose30_dirtyDose(i).doseGrid,dvh_MDirtyDoseOverdosing_DirtyDose30_dirtyDose(i).volumePoints,'Color','black','LineStyle',linestyle(j),'LineWidth',1.5)
    end
end

xlabel('DD in Gy')
xlim = 3;
ylabel('Volume in %')
% legend('Rectum normal','Rectum LETd','Rectum dirty dose','Rectum LETd & Target LETd','Rectum dirtyDose & Target dirtyDose','PTV 68 normal','PTV 68 LETd', ...
%     'PTV 68 dirty dose','PTV 68 LETd & Target LETd','PTV 68 dirtyDose & Target dirtyDose','Bladder normal','Bladder LETd', ...
%     'Bladder dirty dose','Bladder LETd & Target LETd','Bladder dirtyDose & Target dirtyDose','Body normal','Body LETd','Body dirty dose', ...
%     'Body LETd & Target LETd','Body dirtyDose & Target dirtyDose',Location='eastoutside')

% subplot(3,1,3)
% for i = [1 6 8 9]
%     for j = [1 2 3 4]
%         plot(dvh_normalMLETd(i).doseGrid,dvh_normalMLETd(i).volumePoints,'Color','red','LineStyle',linestyle(j),'LineWidth',1.5)
%         hold on
%         plot(dvh_mLETDoseOverdosing_LETd(i).doseGrid,dvh_mLETDoseOverdosing_LETd(i).volumePoints,'Color','blue','LineStyle',linestyle(j),'LineWidth',1.5)
%         hold on
%         plot(dvh_MDirtyDoseOverdosing_LETd(i).doseGrid,dvh_MDirtyDoseOverdosing_LETd(i).volumePoints,'Color','green','LineStyle',linestyle(j),'LineWidth',1.5)
%         hold on
%         plot(dvh_mLETDoseOverdosing_mLETDose120_LETd(i).doseGrid,dvh_mLETDoseOverdosing_mLETDose120_LETd(i).volumePoints,'Color','magenta','LineStyle',linestyle(j),'LineWidth',1.5)
%         hold on
%         plot(dvh_MDirtyDoseOverdosing_DirtyDose30_LETd(i).doseGrid,dvh_MDirtyDoseOverdosing_DirtyDose30_LETd(i).volumePoints,'Color','black','LineStyle',linestyle(j),'LineWidth',1.5)
%     end
% end

% xlim = 10;
% xlabel('LETd in keV/µm')
% ylabel('Volume in %')
% 
% fig = gcf;

% add legend
Lgnd = legend('Rectum Ref','Rectum LETxD__ SO','Rectum DD__ SO','PTV 68 Ref','PTV 68 LETxD__ SO', ...
    'PTV 68 DD__ SO','Bladder Ref','Bladder LETxD__ SO', ...
    'Bladder DD__ SO',Location='southoutside',NumColumns=1);

Lgnd.Position(1) = -0.012;
Lgnd.Position(2) = 0.01;

%%
resultGUIRefNew = resultGUI;
%%
% phys_P10_L4 = matRad_calcQualityIndicators(cst,pln,resultGUI_p10_L4.physicalDose);
RBExD_DD2 = matRad_calcQualityIndicators(cst,pln,resultGUI.RBExD);
% LET_P10_L4 = matRad_calcQualityIndicators(cst,pln,resultGUI_p10_L4.LETd);

%%
figure
subplot(3,3,1)
cube = resultGUI_MDirtyDoseOverdosing.physicalDose;
plane = 3;
slice = 35;
doseWindow = [0 2.5];

matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('(a) dirty dose overdosing (physical dose)')
zoom(1.3)

subplot(3,3,4)
cube = resultGUI_MDirtyDoseOverdosing.dirtyDose;
plane = 3;
slice = 35;
doseWindow = [0 1];

matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('(d) dirty dose overdosing (dirty dose)')
zoom(1.3)

subplot(3,3,7)
cube = resultGUI_MDirtyDoseOverdosing.LETd;
plane = 3;
slice = 35;
doseWindow = [0 3];

matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('(g) dirty dose overdosing (LETd)')
zoom(1.3)

subplot(3,3,2)
cube = resultGUI_mLETDoseOverdosing.physicalDose;
plane = 3;
slice = 35;
doseWindow = [0 2.5];

matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('(b) mLETDose overdosing (physical dose)')
zoom(1.3)

subplot(3,3,5)
cube = resultGUI_mLETDoseOverdosing.dirtyDose;
plane = 3;
slice = 35;
doseWindow = [0 1];

matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('(e) mLETDose overdosing (dirty dose)')
zoom(1.3)

subplot(3,3,8)
cube = resultGUI_mLETDoseOverdosing.LETd;
plane = 3;
slice = 35;
doseWindow = [0 3];

matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('(h) mLETDose overdosing (LETd)')
zoom(1.3)
% 
% cube = resultGUI_M2UClusterLETd.LETd;
% plane = 3;
% slice = 35;
% doseWindow = [0 3];
% 
% figure,
% matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
% title('M2U Cluster LETd (LETd)')
% zoom(1.3)

subplot(3,3,3)
cube = resultGUI_normalM.physicalDose;
plane = 3;
slice = 35;
doseWindow = [0 2.5];

matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('(c) original plan (physical dose)')
zoom(1.3)

subplot(3,3,6)
cube = resultGUI_normalM.dirtyDose;
plane = 3;
slice = 35;
doseWindow = [0 1];

matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('(f) original plan (dirty dose)')
zoom(1.3)

subplot(3,3,9)
cube = resultGUI_normalM.LETd;
plane = 3;
slice = 35;
doseWindow = [0 3];

matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('(i) original plan (LETd)')
zoom(1.3)

cube = resultGUI_withoutM2U.LETd;
plane = 3;
slice = 35;
doseWindow = [0 3];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('M2U mLETDose without cluster (LETd)')
zoom(1.3)

cube = resultGUI_withoutM2U.LETd - resultGUI_M2UClusterDirtyDose.LETd;
plane = 3;
slice = 35;
doseWindow = [-3 3];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference M2U mLETDose and M2U Cluster Dirty Dose')
zoom(1.3)
% 
% cube = resultGUI_p100_L5.RBExD - resultGUI_p70_L5.RBExD;
% plane = 3;
% slice = 80;
% doseWindow = [-0.2 0.2];
% 
% figure,
% matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
% title('difference plan P100 minus P70')
% zoom(1.3)


% maximumRef = max(resultGUIRef.RBExD(cst{9,4}{1}));
% minimumRef = min(resultGUIRef.RBExD(cst{9,4}{1}));

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 100 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cube = zeros(183,183,90);
cube(cst{6,4}{1}) = 1;
vResolution = ct.resolution;
vMargin = [];
vMargin.x = 5;
vMargin.y = 5;
vMargin.z = 5;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);

cst{11,1}    = 10;
cst{11,2}    = 'Margin';
cst{11,3}    = 'OAR';
cst{11,4}{1} = find(mVOIEnlarged);
cst{11,5}    = cst{3,5};

cst{11,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,30)); 

cube = zeros(183,183,90);
cube(cst{4,4}{1}) = 1;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);

cst{12,1} = 11;
cst{12,2} = 'Margin';
cst{12,3} = 'OAR';
cst{12,4}{1} = find(mVOIEnlarged);
cst{12,5} = cst{11,5};

cst{12,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,30));

cube = zeros(183,183,90);
cube(cst{10,4}{1}) = 1;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);

cst{13,1} = 12;
cst{13,2} = 'Margin';
cst{13,3} = 'OAR';
cst{13,4}{1} = find(mVOIEnlarged);
cst{13,5} = cst{11,5};

cst{13,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,30));

cube = zeros(183,183,90);
cube(cst{6,4}{1}) = 1;
cube(cst{1,4}{1}) = 1;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);

cst{14,1} = 13;
cst{14,2} = 'Margin';
cst{14,3} = 'OAR';
cst{14,4}{1} = find(mVOIEnlarged);
cst{14,5} = cst{11,5};

cst{14,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,30));

cst{9,6}{2} = struct(DoseObjectives.matRad_SquaredOverdosing(300,maximumRef));
cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(10,10));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 120 240 270];
pln.propStf.couchAngles   = [0 0 0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultTarget = resultGUI;

cube = resultTarget.RBExD;
plane = 3;
slice = 34;
% doseWindow = [0 max(cube(:))];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
zoom(1.3)

cube = resultTarget.dirtyDose;
plane = 3;
slice = 34;
% doseWindow = [0 max(cube(:))];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
zoom(1.3)

%% Differences

cube = resultGUIRef.RBExD - resultTarget.RBExD;
plane = 3;
slice = 34;
doseWindow = [-1 1];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')
zoom(1.3)

cube = resultGUIRef.dirtyDose - resultTarget.dirtyDose;
plane = 3;
slice = 34;
% doseWindow = [0 max(cube(:))];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')
zoom(1.3)
