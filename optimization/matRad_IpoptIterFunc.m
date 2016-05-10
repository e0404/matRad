function flag = matRad_IpoptIterFunc(iter,objective,~,~,cst,ct)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: iter function
% 
% call
%   Flag = matRad_IpoptIterFunc(iter,objective,parameter,maxiter,figureWait)
%
% input
%   iter:       current number of iteration
%   objective:  current value of objective
%   parameter:  struct with current parameter of optimization
%   maxiter:    max number of iterations
%   figureWait: handle to waitbar
%
% output
%   Flag: indicates if optimization stops (false) / proceeds (true)
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('E:\Mescher\15_DCH_objectiv_tests');

global matRad_STRG_C_Pressed
global matRad_objective_function_value;
global matRad_iteration;
global matRad_global_d;

% update global objective function value
matRad_objective_function_value(iter+1) = objective;

% update global iteration
matRad_iteration = iter + 1;

% check if user pressed ctrl-c
if matRad_STRG_C_Pressed
    flag = false;
else
    flag = true;
end

% plot objective function output
figHandles = get(0,'Children');
if ~isempty(figHandles)
    IdxHandle = strcmp(get(figHandles,'Name'),'Progress of Optimization');
else
    IdxHandle = [];
end

if any(IdxHandle)
    figOpt = figHandles(IdxHandle);
    %AxesInfigOpt = findall(figOpt,'type','axes');
    AxesInfigOpt = subplot(1,2,1);
    set(AxesInfigOpt,'NextPlot', 'replacechildren')
    children = get(AxesInfigOpt,'children');
    delete(children);
else
    figOpt = figure('Name','Progress of Optimization','NumberTitle','off','Color',[.5 .5 .5]);
    %hold on, grid on, grid minor,
    %AxesInfigOpt = findall(figOpt,'type','axes');
    AxesInfigOpt = subplot(1,2,1);
    hold on, grid on, grid minor,
end

defaultFontSize = 14;
set(AxesInfigOpt,'YScale','log');
title(AxesInfigOpt,'Progress of Optimization','LineWidth',defaultFontSize),
xlabel(AxesInfigOpt,'# iterations','Fontsize',defaultFontSize),ylabel(AxesInfigOpt,'objective function value','Fontsize',defaultFontSize)

% draw updated axes
plot(AxesInfigOpt,0:1:iter,matRad_objective_function_value,'xb','LineWidth',1.5);
subplot(1,2,2)
plot_CT_Dose_VOIs(ct,cst([6 11],:),reshape(matRad_global_d{1},ct.cubeDim),2,90,0,1,0,0)
drawnow

end