function Flag = matRad_IpoptIterFunc(iter,objective,parameter,maxiter,figureWait)
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

global matRad_STRG_C_Pressed
global matRad_objective_function_value;

% update global objective function value
matRad_objective_function_value(iter+1) = objective;

% update waitbar
waitbar(iter/maxiter,figureWait,['Optimizing beam intensities (iter = ',num2str(iter),')']);

% check if user pushed cancel button
if getappdata(figureWait,'canceling')
    Flag = false;
else
    Flag = true;
end

% check if user pressed ctrl-c
if matRad_STRG_C_Pressed
    Flag = false;
end

% plot objective function output
try
    figHandles = get(0,'Children');
    IdxHandle = [];
    if ~isempty(figHandles)
        v=version;
        if str2num(v(1:3))>=8.5
            IdxHandle = strcmp({figHandles(:).Name},'Progress of Optimization');
        else
            IdxHandle = strcmp(get(figHandles,'Name'),'Progress of Optimization');
        end
    end
    if  any(IdxHandle)
        figOpt = figHandles(IdxHandle);
        AxesInfigOpt = findall(figOpt,'type','axes');
        set(AxesInfigOpt,'NextPlot', 'replacechildren')
        v=version;
        if str2num(v(1:3))>=8.5
            delete(AxesInfigOpt.Children);
        else
            children = get(AxesInfigOpt,'children');
            delete(children);
        end
    else
        figOpt=figure('Name','Progress of Optimization','NumberTitle','off');
        hold on, grid on, grid minor,
        AxesInfigOpt = findall(figOpt,'type','axes');
    end
    % prevent closure of window and show busy state
    set(figOpt,'CloseRequestFcn','');
    set(figOpt,'pointer','watch');

    set(AxesInfigOpt,'YScale','log');
    title(AxesInfigOpt,'Progress of Optimization','LineWidth',14),
    xlabel(AxesInfigOpt,'# iterations','Fontsize',14),ylabel(AxesInfigOpt,'objective function value','Fontsize',14)
catch 
    warning('couldnt initialize figure to plot the objective value')
end

% draw updated axes
axes(AxesInfigOpt)
plot(AxesInfigOpt,0:1:iter,matRad_objective_function_value,'xb','LineWidth',1.5);
drawnow

% revert busy state and enable close button
set(figOpt,'CloseRequestFcn',get(0,'DefaultFigureCloseRequestFcn'));
set(figOpt,'pointer','arrow');


end