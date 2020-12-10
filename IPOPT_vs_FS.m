load TG119.mat

pln.radiationMode = 'photons';  
pln.machine       = 'Generic';
pln.propOpt.bioOptimization = 'none';    
pln.numOfFractions         = 1;
pln.propStf.gantryAngles   = [0:360/7:359];
%pln.propStf.gantryAngles   = [0];
pln.propStf.couchAngles    = zeros(1, numel(pln.propStf.gantryAngles));
pln.propStf.bixelWidth     = 5;
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
%pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
%pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
%pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution = ct.resolution;
pln.propOpt.runSequencing = 1;
pln.propOpt.runDAO        = 0;
stf                      = matRad_generateStf(ct,cst,pln);
dij = matRad_calcPhotonDose(ct,stf,pln,cst);
%pln.propOpt.tol_obj = 1e-6;
%pln.propOpt.tol_violation = 1e-6;
%pln.propOpt.accepted_violation = 1e-5;

cst{1, 5}.visibleColor = [0.5 0.5 0.5];
cst{2, 5}.visibleColor = [0 0 0];
cst{3, 5}.visibleColor = [0.4 0.4470 0.7410];


%Plot settings
plane      = 3;
doseWindow = [0 70];
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);


%% Set Optimization

cst(:,6) = [];

%Core
cst{1, 6}{1}=DoseObjectives.matRad_SquaredOverdosing(100, 20)

%Target
cst{2, 6}{1}=DoseObjectives.matRad_SquaredDeviation(1000, 60);
cst{2, 6}{2}=DoseConstraints.matRad_MinMaxDose(59, 61, 1, 1.9);

%Body
cst{3, 6}{1}=DoseObjectives.matRad_SquaredOverdosing(30, 30);
%cst{1, 6}{1}=DoseConstraints.matRad_MinMaxDose(0, 20, 1, 0.9);


%% IPOPT
pln.propOpt.optimizer = "IPOPT";
%pln.propOpt.feasibility_seeker = "AMS_sim";
pln.propOpt.max_iter = 750;
pln.propOpt.max_time = 3600;
%pln.propOpt.lambda = 1;
%pln.propOpt.weighted = true;
tic;
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
time = toc;


%% Plot IPOPT

hfPlan = figure; 
subplot(2,2,1);
title('IPOPT plan');
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],0.75,[],[],doseWindow,[],[],[],[],'LineWidth',2);

xlim([30 138]);
ylim([55 115]);

usedOpt = resultGUI.usedOptimizer;

hFvals = figure;
title('Convergence');
subplot(2,2,1);plot(0:numel(usedOpt.allObjectiveFunctionValues)-1,usedOpt.allObjectiveFunctionValues,'x'); xlabel('# Iteration'); ylabel('Obj. Function'); grid('minor'); set(gca,'YScale','log'); hold on;
subplot(2,2,2);plot(usedOpt.timeIter,usedOpt.allObjectiveFunctionValues,'x'); xlabel('Time [s]'); ylabel('Obj. Function'); grid('minor'); set(gca,'YScale','log'); hold on;
subplot(2,2,3);plot(0:numel(usedOpt.allConstraintViolations)-1,usedOpt.allConstraintViolations,'x'); xlabel('# Iteration'); ylabel('Constr. Violation'); grid('minor'); hold on;
subplot(2,2,4);plot(usedOpt.timeIter,usedOpt.allConstraintViolations,'x'); xlabel('Time [s]'); ylabel('Constr. Violation'); grid('minor'); hold on;



%% Super AMS_sim
pln.propOpt.optimizer = "Superization";
pln.propOpt.feasibility_seeker = "AMS_sequential";
pln.propOpt.max_iter = 750;
pln.propOpt.max_time = 3600;
pln.propOpt.lambda = 1;
pln.propOpt.weighted = true;



tic;
resultGUI_super = matRad_fluenceOptimization(dij,cst,pln);
time = toc;

%% Plot AMS

figure(hfPlan);
subplot(2,2,2);
title('Superiorized plan');
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_super.physicalDose,plane,slice,[],0.75,[],[],doseWindow,[],[],[],[],'LineWidth',2);
xlim([30 138]);
ylim([55 115]);

usedOpt = resultGUI_super.usedOptimizer;

figure(hFvals);
legEntries = {'IPOPT','Superiorization'};
subplot(2,2,1);plot(0:numel(usedOpt.allObjectiveFunctionValues)-1,usedOpt.allObjectiveFunctionValues,'x'); legend(legEntries); 
subplot(2,2,2);plot(usedOpt.timeIter,usedOpt.allObjectiveFunctionValues,'x');legend(legEntries);
subplot(2,2,3);plot(0:numel(usedOpt.allConstraintViolations)-1,usedOpt.allConstraintViolations,'x'); legend(legEntries);
subplot(2,2,4);plot(usedOpt.timeIter,usedOpt.allConstraintViolations,'x'); legend(legEntries); 

%% 


absDiffCube = resultGUI.physicalDose-resultGUI_super.physicalDose;

diffMax = max(abs(absDiffCube(:)));
diffWindow = [-diffMax diffMax];

figure(hfPlan);
subplot(2,2,3);
title('IPOPT plan - Sup. Plan');
matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],0.75,[],diffMap,diffWindow,[],[],[],[],'LineWidth',2);
xlim([30 138]);
ylim([55 115]);

%% Obtain dose statistics
% Two more columns will be added to the cst structure depicting the DVH and
% standard dose statistics such as D95,D98, mean dose, max dose etc.
[dvh,qi]               = matRad_indicatorWrapper(cst,pln,resultGUI);
[dvh_super,qi_super] = matRad_indicatorWrapper(cst,pln,resultGUI_super);

%% Show DVHs
figure(hfPlan);
subplot(2,2,4);
title('DVHs');
matRad_showDVH(dvh,cst,pln,1)
hold on
matRad_showDVH(dvh_super,cst,pln,2)

names1 = cellfun(@(c) sprintf('%s - IPOPT',c),cst(:,2),'UniformOutput',false);
names2 = cellfun(@(c) sprintf('%s - Sup.',c),cst(:,2),'UniformOutput',false);

allNames = [names1 names2];

legend(allNames{:},'Location','NorthOutside','NumColumns',2);

