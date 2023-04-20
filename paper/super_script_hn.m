load HEAD_AND_NECK_super.mat

%%Change method

meth = 1;

for i = 1:size(cst,1)
    for o = 1:numel(cst{i,6})
        fun = cst{i,6}{o};
        if strcmp(fun.className,'DoseConstraints.matRad_MinMaxDose')
            cst{i,6}{o}.parameters 
        end
    end
end

%% IPOPT
opti = matRad_OptimizerIPOPT;
opti.options.max_iter = 2000;
opti.options.max_cpu_time = 3600;
opti.options.dual_inf_tol              = 1e2; % (Opt2)
opti.options.constr_viol_tol           = 1e-2; % (Opt3)
opti.options.acceptable_iter           = 5;    % (Acc1)
opti.options.acceptable_tol            = 1e10; % (Acc2) %The scale of our objective function is variable, so we relax this tolerance
opti.options.acceptable_obj_change_tol = 1e-4; % (Acc6), Solved To Acceptable Level if (Acc1),...,(Acc6) fullfiled

pln.propOpt.optimizer = opti;

tic;
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
time = toc;

%% Super AMS_sim
opti = matRad_OptimizerSuperization;
opti.feasibility_seeker = "AMS_sequential";
opti.max_iter = 2000;
opti.max_time = 3600;
opti.lambda = 1;
opti.alpha = 0.99;
opti.weighted = true;
opti.control_sequence = 'weight';
opti.weight_decay = 1;
opti.warm_start = true;
opti.accepted_tol_change = 1e-4;
opti.accepted_violation = 1e-6;
opti.accepted_max_violation = 1e-2;
pln.propOpt.optimizer = opti;

tic;
resultGUI_super = matRad_fluenceOptimization(dij,cst,pln);
time = toc;


%% Plot IPOPT
plane = 3;
slice = 28;
doseWindow = [0 75]./pln.numOfFractions;
xWindow = [30 128];
yWindow = [30 128];


hfPlan = figure; 
hfPlan.WindowState = 'Maximized';
axIpoptPlan = subplot(2,2,1);
title('IPOPT plan');
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],0.75,[],[],doseWindow,[],[],[],[],'LineWidth',2);
xlim(xWindow);
ylim(yWindow);

usedOpt = resultGUI.usedOptimizer;

hFvals = figure;
hFvals.WindowState = 'Maximized';
%title('Convergence');
subplot(2,2,1);plot(0:numel(usedOpt.allObjectiveFunctionValues)-1,usedOpt.allObjectiveFunctionValues,'x'); xlabel('# Iteration'); ylabel('Obj. Function'); grid('minor'); set(gca,'YScale','log'); hold on;
subplot(2,2,2);plot(usedOpt.timeIter,usedOpt.allObjectiveFunctionValues,'x'); xlabel('Time [s]'); ylabel('Obj. Function'); grid('minor'); set(gca,'YScale','log'); hold on;
subplot(2,2,3);plot(0:numel(usedOpt.allConstraintViolations)-1,usedOpt.allConstraintViolations,'x'); xlabel('# Iteration'); ylabel('Constr. Violation'); grid('minor'); hold on;
subplot(2,2,4);plot(usedOpt.timeIter,usedOpt.allConstraintViolations,'x'); xlabel('Time [s]'); ylabel('Constr. Violation'); grid('minor'); hold on;

%% Plot AMS

figure(hfPlan);

axSuperPlan = subplot(2,2,2);
title('Superiorized plan');
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_super.physicalDose,plane,slice,[],0.75,[],[],doseWindow,[],[],[],[],'LineWidth',2);
xlim(xWindow);
ylim(yWindow);


usedOpt = resultGUI_super.usedOptimizer;

figure(hFvals);
legEntries = {'IPOPT','Superiorization'};
axObjVsIter  = subplot(2,2,1); plot(0:numel(usedOpt.allObjectiveFunctionValues)-1,usedOpt.allObjectiveFunctionValues,'x');  legend(legEntries); 
axObjVsTime  = subplot(2,2,2); plot(usedOpt.timeIter,usedOpt.allObjectiveFunctionValues,'x');                               legend(legEntries);
axViolVsIter = subplot(2,2,3) ;plot(0:numel(usedOpt.allConstraintViolations)-1,usedOpt.allConstraintViolations,'x');        legend(legEntries);
axViolVsTime = subplot(2,2,4); plot(usedOpt.timeIter,usedOpt.allConstraintViolations,'x');                                  legend(legEntries); 

%% 


absDiffCube = resultGUI.physicalDose-resultGUI_super.physicalDose;

diffMax = max(abs(absDiffCube(:)));
diffWindow = [-diffMax diffMax];

figure(hfPlan);
axDiffPlan = subplot(2,2,3);
title('IPOPT plan - Sup. Plan');
matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],0.75,[],diffMap,diffWindow,[],[],[],[],'LineWidth',2);
xlim(xWindow);
ylim(yWindow);



%% Obtain dose statistics
% Two more columns will be added to the cst structure depicting the DVH and
% standard dose statistics such as D95,D98, mean dose, max dose etc.
[dvh,qi]               = matRad_indicatorWrapper(cst,pln,resultGUI);
[dvh_super,qi_super] = matRad_indicatorWrapper(cst,pln,resultGUI_super);

%% Show DVHs
figure(hfPlan);
axDVHs = subplot(2,2,4);
title('DVHs');
matRad_showDVH(dvh,cst,pln,1)
hold on
matRad_showDVH(dvh_super,cst,pln,2)

visibleIx = cellfun(@(c) c.Visible == 1,cst(:,5));
names1 = cellfun(@(c) sprintf('%s - IPOPT',c),cst(visibleIx,2),'UniformOutput',false);
names2 = cellfun(@(c) sprintf('%s - Sup.',c),cst(visibleIx,2),'UniformOutput',false);

allNames = [names1 names2];

legend(allNames{:},'Location','NorthOutside','NumColumns',4);

%%
saveAsPngAndFig(hfPlan,'super_hn_plan');
saveAsPngAndFig(hfPlan,'super_hn_vals');


saveAsPngAndFig(axObjVsIter,'super_hn_objVsIter');
saveAsPngAndFig(axObjVsTime,'super_hn_objVsTime');
saveAsPngAndFig(axViolVsIter,'super_hn_violVsIter');
saveAsPngAndFig(axViolVsTime,'super_hn_violVsTime');

saveAsPngAndFig(axSuperPlan,'super_hn_superPlan');
saveAsPngAndFig(axIpoptPlan,'super_hn_ipoptPlan');
saveAsPngAndFig(axDiffPlan,'super_hn_diffPlan');

saveAsPngAndFig(axDVHs,'super_hn_dvhs');

%%
clear dij;
save('super_hn.mat');

%% matlab2tikz
plane = 3;
slice = 28;
doseWindow = [0 75]./pln.numOfFractions;
xWindow = [40 118];
yWindow = [30 108];


hfPlan = figure; 
hfPlan.WindowState = 'Maximized';
%axIpoptPlan = subplot(2,2,1);
set(hfPlan,'DefaultLegendInterpreter','none');
[hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],0.75,[],[],doseWindow,[],[],[],true,'LineWidth',2);
xlim(xWindow);
ylim(yWindow);
title('');
xlabel('')
ylabel('')
xticks([])
yticks([])
set(hCMap.Label,'String','dose [Gy]')
cleanfigure; matlab2tikz('D:\DKFZHomeOffice\paper\superiorization\tikz\super_script_hn_ipoptPlan.tikz','showInfo',true,'width','\figurewidth','height','\figureheight','relativeDataPath','tikz/');
close(hfPlan);

%%
hfPlan = figure; 
hfPlan.WindowState = 'Maximized';
%axIpoptPlan = subplot(2,2,1);
set(hfPlan,'DefaultLegendInterpreter','none');
[hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_super.physicalDose,plane,slice,[],0.75,[],[],doseWindow,[],[],[],true,'LineWidth',2);
xlim(xWindow);
ylim(yWindow);
title('');
xlabel('')
ylabel('')
xticks([])
yticks([])
set(hCMap.Label,'String','dose [Gy]')
cleanfigure; matlab2tikz('D:\DKFZHomeOffice\paper\superiorization\tikz\super_script_hn_superPlan.tikz','showInfo',true,'width','\figurewidth','height','\figureheight','relativeDataPath','tikz/');
close(hfPlan);

%%
hfPlan = figure; 
hfPlan.WindowState = 'Maximized';
%axIpoptPlan = subplot(2,2,1);
set(hfPlan,'DefaultLegendInterpreter','none');
[hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],0.75,[],diffMap,diffWindow,[],[],[],false,'LineWidth',2);
xlim(xWindow);
ylim(yWindow);
title('');
xlabel('')
ylabel('')
xticks([])
yticks([])
set(hCMap.Label,'String','dose [Gy]')
cleanfigure; matlab2tikz('D:\DKFZHomeOffice\paper\superiorization\tikz\super_script_hn_diffPlan.tikz','showInfo',true,'width','\figurewidth','height','\figureheight','relativeDataPath','tikz/');
close(hfPlan);

%%
hfPlan = figure;
%axDVHs = subplot(2,2,4);
%title('DVHs');
notVisIx = [4 6 8 9 10 11 12 20 21 22 23 18 1];
for i = 1:numel(notVisIx)
    cst{notVisIx(i),5}.Visible = false;
end
matRad_showDVH(dvh,cst,pln,1)
hold on
matRad_showDVH(dvh_super,cst,pln,2)

visibleIx = cellfun(@(c) c.Visible == 1,cst(:,5));
names1 = cellfun(@(c) sprintf('%s - IPOPT',c),cst(visibleIx,2),'UniformOutput',false);
names2 = cellfun(@(c) sprintf('%s - Sup.',c),cst(visibleIx,2),'UniformOutput',false);

allNames = [names1 names2];
xlim([0 80]);

legend(allNames{:},'Location','NorthOutside','NumColumns',4);
cleanfigure; matlab2tikz('D:\DKFZHomeOffice\paper\superiorization\tikz\super_script_hn_dvhs.tikz','showInfo',true,'width','\figurewidth','height','\figureheight');

%%
hfPlan = figure;
usedOpt = resultGUI.usedOptimizer;
semilogy(0:numel(usedOpt.allObjectiveFunctionValues)-1,usedOpt.allObjectiveFunctionValues,'x'); hold on;
usedOpt = resultGUI_super.usedOptimizer;
semilogy(0:numel(usedOpt.allObjectiveFunctionValues)-1,usedOpt.allObjectiveFunctionValues,'x'); 
grid minor;
xlabel('#Iterations');
ylabel('Objective function value');
legend({'Optimization','Superiorization'});
cleanfigure; matlab2tikz('D:\DKFZHomeOffice\paper\superiorization\tikz\super_script_hn_objvsiter.tikz','showInfo',true,'width','\figurewidth','height','\figureheight');

%%
hfPlan = figure;
usedOpt = resultGUI.usedOptimizer;
semilogy(usedOpt.timeIter,usedOpt.allObjectiveFunctionValues,'x'); hold on;
usedOpt = resultGUI_super.usedOptimizer;
semilogy(usedOpt.timeIter,usedOpt.allObjectiveFunctionValues,'x'); 
grid minor;
xlabel('time [s]');
ylabel('Objective function value');
legend({'Optimization','Superiorization'})
cleanfigure; matlab2tikz('D:\DKFZHomeOffice\paper\superiorization\tikz\super_script_hn_objvstime.tikz','showInfo',true,'width','\figurewidth','height','\figureheight');

%%
hfPlan = figure;
usedOpt = resultGUI.usedOptimizer;
plot(0:numel(usedOpt.allConstraintViolations)-1,usedOpt.allConstraintViolations,'x'); hold on;
usedOpt = resultGUI_super.usedOptimizer;
plot(0:numel(usedOpt.allConstraintViolations)-1,usedOpt.allConstraintViolations,'x'); 
grid minor;
xlabel('#Iterations');
ylabel('Constraint violation');
legend({'Optimization','Superiorization'});
cleanfigure; matlab2tikz('D:\DKFZHomeOffice\paper\superiorization\tikz\super_script_hn_constrvsiter.tikz','showInfo',true,'width','\figurewidth','height','\figureheight');

%%
hfPlan = figure;
usedOpt = resultGUI.usedOptimizer;
plot(usedOpt.timeIter,usedOpt.allConstraintViolations,'x'); hold on;
usedOpt = resultGUI_super.usedOptimizer;
plot(usedOpt.timeIter,usedOpt.allConstraintViolations,'x'); 
grid minor;
xlabel('time [s]');
ylabel('Constraint violation');
legend({'Optimization','Superiorization'})
cleanfigure; matlab2tikz('D:\DKFZHomeOffice\paper\superiorization\tikz\super_script_hn_constrvstime.tikz','showInfo',true,'width','\figurewidth','height','\figureheight');
