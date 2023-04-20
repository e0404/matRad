clear
load TG119_super.mat
%load TG119_super_protons.mat
matRad_rc
%% Set Optimization

cst(:,6) = [];

%Core
cst{1, 6}{1}=DoseObjectives.matRad_SquaredOverdosing(100, 20)
%cst{1, 6}{2}=DoseConstraints.matRad_MinMaxDose(0, 30, 1, 0.9);

%Target
cst{2, 6}{1}=DoseObjectives.matRad_SquaredOverdosing(1000, 61);
cst{2, 6}{2}=DoseObjectives.matRad_SquaredUnderdosing(1000, 59);
%cst{2, 6}{2}=DoseConstraints.matRad_MinMaxDose(59, 61, 1, 1.9);

%Body
cst{3, 6}{1}=DoseObjectives.matRad_SquaredOverdosing(30, 30);
%cst{1, 6}{1}=DoseConstraints.matRad_MinMaxDose(0, 20, 1, 0.9);


%% IPOPT
opti = matRad_OptimizerIPOPT;
opti.options.max_iter = 1000;
opti.options.max_cpu_time = 3600;
opti.options.dual_inf_tol              = 1e-2; % (Opt2)
opti.options.constr_viol_tol           = 1e-2; % (Opt3)
opti.options.acceptable_iter           = 5;    % (Acc1)
opti.options.acceptable_tol            = 1e10; % (Acc2) %The scale of our objective function is variable, so we relax this tolerance
opti.options.acceptable_obj_change_tol = 1e-4; % (Acc6), Solved To Acceptable Level if (Acc1),...,(Acc6) fullfiled

pln.propOpt.optimizer = opti;

tic;
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
time = toc;


%% Set Optimization

cst(:,6) = [];

%Core
cst{1, 6}{1}=DoseConstraints.matRad_MinMaxDose(0, 20, 2, 100);
%cst{1, 6}{2}=DoseObjectives.matRad_SquaredOverdosing(100, 20);

%Target
cst{2, 6}{1}=DoseConstraints.matRad_MinMaxDose(59, 61, 2, 1000);

%cst{2, 6}{2}=DoseObjectives.matRad_SquaredOverdosing(1000, 61);
%cst{2, 6}{3}=DoseObjectives.matRad_SquaredUnderdosing(1000, 59);

%Body
cst{3, 6}{1}=DoseConstraints.matRad_MinMaxDose(0, 30, 2, 30);
%cst{3, 6}{2}=DoseObjectives.matRad_SquaredOverdosing(30, 30);

%cst{3, 6}{1}=DoseConstraints.matRad_MinMaxDose(0, 30, 1, 50);


%% Super AMS_sim
opti = matRad_OptimizerSuperization;
opti.feasibility_seeker = "AMS_sequential";
opti.max_iter = 1000;
opti.max_time = 3600;
opti.lambda = 1.5;
opti.weighted = true;
opti.control_sequence = 'weight';
opti.weight_decay = 0.99; %0.99;
opti.warm_start = true;
opti.ignoreObjective = true;
pln.propOpt.optimizer = opti;

tic;
resultGUI_super = matRad_fluenceOptimization(dij,cst,pln);
time = toc;


%% Plot IPOPT

doseWindow = [0 75];
%isoLevels = [0.2 0.4 0.6 0.8 0.95 1.0 1.05 1.1]*60;
isoLevels = [];

hfPlan = figure('WindowState','maximized'); 
axIpoptPlan = subplot(2,2,1);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],0.75,[],[],doseWindow,isoLevels,[],[],[],'LineWidth',2);
%title('IPOPT plan');

xlim([30 138]);
ylim([55 115]);

usedOpt = resultGUI.usedOptimizer;

hFvals = figure('WindowState','maximized');
subplot(2,2,1);plot(0:numel(usedOpt.allObjectiveFunctionValues)-1,usedOpt.allObjectiveFunctionValues,'x'); xlabel('# Iteration'); ylabel('Obj. Function'); grid('minor'); set(gca,'YScale','log'); hold on;
subplot(2,2,2);plot(usedOpt.timeIter,usedOpt.allObjectiveFunctionValues,'x'); xlabel('Time [s]'); ylabel('Obj. Function'); grid('minor'); set(gca,'YScale','log'); hold on;
subplot(2,2,3);plot(0:numel(usedOpt.allConstraintViolations)-1,usedOpt.allConstraintViolations,'x'); xlabel('# Iteration'); ylabel('Constr. Violation'); grid('minor'); hold on;
subplot(2,2,4);plot(usedOpt.timeIter,usedOpt.allConstraintViolations,'x'); xlabel('Time [s]'); ylabel('Constr. Violation'); grid('minor'); hold on;
title('Convergence');

%% Plot AMS

figure(hfPlan);
axSuperPlan = subplot(2,2,2);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_super.physicalDose,plane,slice,[],0.75,[],[],doseWindow,isoLevels,[],[],[],'LineWidth',2);
xlim([30 138]);
ylim([55 115]);

%title('AMS Plan');

usedOpt = resultGUI_super.usedOptimizer;

figure(hFvals);
legEntries = {'IPOPT','Superiorization'};
axObjVsIter  = subplot(2,2,1);plot(0:numel(usedOpt.allObjectiveFunctionValues)-1,usedOpt.allObjectiveFunctionValues,'x'); legend(legEntries); 
axObjVsTime  = subplot(2,2,2);plot(usedOpt.timeIter,usedOpt.allObjectiveFunctionValues,'x');legend(legEntries);
axViolVsIter = subplot(2,2,3);plot(0:numel(usedOpt.allConstraintViolations)-1,usedOpt.allConstraintViolations,'x'); legend(legEntries);
axViolVsTime = subplot(2,2,4);plot(usedOpt.timeIter,usedOpt.allConstraintViolations,'x'); legend(legEntries); 

%% 


absDiffCube = resultGUI.physicalDose-resultGUI_super.physicalDose;

diffMax = max(abs(absDiffCube(:)));
diffWindow = [-diffMax diffMax];

figure(hfPlan);
axDiffPlan = subplot(2,2,3);
%title('IPOPT plan - AMS Plan');
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
axDVHs = subplot(2,2,4);
%title('DVHs');
matRad_showDVH(dvh,cst,pln,1)
hold on
matRad_showDVH(dvh_super,cst,pln,2)
xlim(doseWindow);

visibleIx = cellfun(@(c) c.Visible == 1,cst(:,5));
names1 = cellfun(@(c) sprintf('%s - IPOPT',c),cst(visibleIx,2),'UniformOutput',false);
names2 = cellfun(@(c) sprintf('%s - Sup.',c),cst(visibleIx,2),'UniformOutput',false);

allNames = [names1 names2];

legend(allNames{:},'Location','NorthOutside','NumColumns',2);

%% Save
saveAsPngAndFig(hfPlan,'super_script_basic_plan');
saveAsPngAndFig(hFvals,'super_script_basic_vals');


saveAsPngAndFig(axObjVsIter,'super_script_basic_objVsIter');
saveAsPngAndFig(axObjVsTime,'super_script_basic_objVsTime');
saveAsPngAndFig(axViolVsIter,'super_script_basic_violVsIter');
saveAsPngAndFig(axViolVsTime,'super_script_basic_violVsTime');

saveAsPngAndFig(axSuperPlan,'super_script_basic_superPlan');
saveAsPngAndFig(axIpoptPlan,'super_script_basic_ipoptPlan');
saveAsPngAndFig(axDiffPlan,'super_script_basic_diffPlan');

saveAsPngAndFig(axDVHs,'supe_script_basic_dvhs')

%%
close all;
save('super_script_basic.mat');