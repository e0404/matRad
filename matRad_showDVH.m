function matRad_showDVH(cst,pln,scenIx,lineStyleIndicator)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad dvh calculation
% 
% call
%   matRad_calcDVH(d,cst,lineStyleIndicator)
%
% input
%   result:             result struct from fluence optimization/sequencing
%   cst:                matRad cst struct
%   scenIx: (optional)  index of scenario (default: 1)
%   lineStyleIndicator: integer (1,2,3,4) to indicate the current linestyle
%                       (hint: use different lineStyles to overlay
%                       different dvhs)
%
% output
%   graphical display of DVH & dose statistics in console   
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

% create new figure and set default line style indicator if not explictly
% specified
if ~exist('scenIx','var') || isempty(scenIx)
    scenIx = 1;
end

if ~exist('lineStyleIndicator','var') || isempty(lineStyleIndicator)
    f = figure('Name','DVH','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
    hold on
    lineStyleIndicator = 1;
else
    hold on
end

numOfVois = size(cst,1);

% Create the column and row names in cell arrays 
cnames = {'dummy_a'};
rnames = cst(:,2);
% Create the uitable
table = uitable(gcf,'Data',zeros(length(rnames),length(cnames)),...
            'ColumnName',cnames,... 
            'RowName',rnames,'ColumnWidth',{70});
        
%% print the dvh
colorMx    = colorcube;
colorMx    = colorMx(1:floor(64/numOfVois):64,:);

lineStyles = {'-',':','--','-.'};

maxDVH = 0;

for i = 1:numOfVois
    if cst{i,5}.Visible
        dvh = cst{i,8}{scenIx};
        % cut off at the first zero value where there is no more signal
        % behind
        ix = find(dvh(2,:)>0,1,'last');
        dvh = dvh(:,1:ix);
        subplot(211);
        plot(dvh(1,:),dvh(2,:),'LineWidth',4,'Color',colorMx(i,:), ...
            'LineStyle',lineStyles{lineStyleIndicator},'DisplayName',cst{i,2});hold on
        maxDVH = max(maxDVH,max(dvh(2,:)));
    end
end

fontSizeValue = 14;
myLegend = legend('show','location','NorthEast');
set(myLegend,'FontSize',10,'Interpreter','none');
legend boxoff


ylim([0 1.1*maxDVH]);
xlim([0 1.2*max(dvh(1,:))]);

grid on,grid minor
box(gca,'on');
set(gca,'LineWidth',1.5,'FontSize',fontSizeValue);
ylabel('Volume [%]','FontSize',fontSizeValue)

if strcmp(pln.bioOptimization,'none')
     xlabel('Dose [Gy]','FontSize',fontSizeValue);
else
     xlabel('RBE x Dose [Gy(RBE)]','FontSize',fontSizeValue);
end

pos = get(subplot(2,1,2),'position');
ylabel('VOIs');
xlabel('dose statistics');
set(subplot(2,1,2),'yTick',[])
set(subplot(2,1,2),'xTick',[])

set(table,'units','normalized')
set(table,'position',pos)

% get quality indicators and fill table
QI = [];
for i = 1:numOfVois
    QI = [QI; cst{i,9}{scenIx}];
end

% remove underscore from display
indicatorNames = fieldnames(QI);
for i = 1:numel(indicatorNames)
    ix = find(indicatorNames{i}(4:end) == '_');
    if ~isempty(ix)
        indicatorNames{i}(ix+3) = '.';
    end
end

set(table,'ColumnName',indicatorNames);
set(table,'Data',(squeeze(struct2cell(QI)))');

