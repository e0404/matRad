function flag = matRad_IpoptIterFunc(iter,objective,~,~)
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
global matRad_iteration;

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
    AxesInfigOpt = findall(figOpt,'type','axes');
    set(AxesInfigOpt,'NextPlot', 'replacechildren')
    children = get(AxesInfigOpt,'children');
    
    if length(children) > 1
        for i = 1:length(children)
            delete(children{i});
        end
        AxesInfigOpt = AxesInfigOpt(end);
    else
        delete(children);
    end

else
    figOpt = figure('Name','Progress of Optimization','NumberTitle','off','Color',[.5 .5 .5]);
    %subplot(2,3,1)
    hold on, grid on, grid minor,
    AxesInfigOpt = findall(figOpt,'type','axes');
end
% ensure to bring optimization window to front also for a re-optimization 
if isdeployed
    figure(figOpt);
end 
defaultFontSize = 14;
set(AxesInfigOpt,'YScale','log');
title(AxesInfigOpt,'Progress of Optimization','LineWidth',defaultFontSize),
xlabel(AxesInfigOpt,'# iterations','Fontsize',defaultFontSize),ylabel(AxesInfigOpt,'objective function value','Fontsize',defaultFontSize)

% draw updated axes
plot(AxesInfigOpt,0:1:iter,matRad_objective_function_value,'xb','LineWidth',1.5);

% draw additional variables
% DVH and DCH scaling parameter
% global kDVH
% global kDCH
% if ~isempty(kDVH) & ~isempty(kDCH)
%     colors = {'b','r','k','g'};
%     for i = 1:size(kDCH,1)
%         h1 = subplot(2,3,2);
%         hold on, grid on, grid minor
%         plot(0:1:iter,kDVH(i,:),'x','Color',colors{i},'LineWidth',1.5)
%         set(h1,'YScale','log');
%         title('kDVH')
% 
%         h2 = subplot(2,3,3);
%         hold on, grid on, grid minor
%         plot(0:1:iter,kDCH(i,:),'x','Color',colors{i},'LineWidth',1.5)
%         set(h2,'YScale','log');
%         title('kDCH')
% 
%     end
% end
% 
% % unscaled constraint value
% global CONSTRAINT
% if ~isempty(CONSTRAINT)
%     colors = {'b','r','k','g'};
%     for i = 1:size(CONSTRAINT,1)
%         h3 = subplot(2,3,4);
%         hold on, grid on, grid minor
%         plot(0:1:iter,CONSTRAINT(i,:),'x','Color',colors{i},'LineWidth',1.5)
%     end
% set(h3,'YScale','lin');
% title('unscaled constraint')
% end
% 
% % min/max jacobian
% colors = {'b','r','k','g'};
% marker = {'v','^','x','*'};
% global JACOBIAN
% if size(JACOBIAN,3) == matRad_iteration & ~isempty(JACOBIAN)
%     
%     for i = 1:min(size(JACOBIAN,1),2)
%         h1 = subplot(2,3,5);
%         hold on, grid on, grid minor
%         %plot(1:size(JACOBIAN,2),abs(JACOBIAN(i,:,matRad_iteration)),'Color',colors{i},'LineWidth',1.5)
%         plot(0:1:iter,squeeze(JACOBIAN(i,1,:)),marker{i},'Color',colors{i},'LineWidth',1.5)
%         plot(0:1:iter,squeeze(JACOBIAN(i,2,:)),marker{i},'Color',colors{i},'LineWidth',1.5)
%         set(h1,'YScale','log');
%         title('minmax(abs(jacobian))')
%         ylim([1e-10 1e10])
%         
%     end
% end
% 
% % min/max gradient
% global GRADIENT
% if size(GRADIENT,3) == matRad_iteration & ~isempty(GRADIENT)
%     
%         h2 = subplot(2,3,6);
%         hold on, grid on, grid minor
%         %plot(1:size(GRADIENT,2),abs(GRADIENT(:,:,matRad_iteration)),'Color',colors{i},'LineWidth',1.5)
%         plot(0:1:iter,squeeze(GRADIENT(1,1,:)),'x','Color',colors{1},'LineWidth',1.5)
%         plot(0:1:iter,squeeze(GRADIENT(1,2,:)),'x','Color',colors{1},'LineWidth',1.5)
%         set(h2,'YScale','log');
%         title('minmax(abs(gradient))')
%         ylim([1e-10 1e10])
% end
% drawnow


end