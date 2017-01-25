function matRad_calcDVH(result,cst,pln,lineStyleIndicator)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad dvh calculation
% 
% call
%   matRad_calcDVH(d,cst,lineStyleIndicator)
%
% input
%   result:             result struct from fluence optimization/sequencing
%   cst:                matRad cst struct
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
if nargin < 4 
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
        
%% calculate and print the dvh
colorMx    = colorcube;
colorMx    = colorMx(1:floor(64/numOfVois):64,:);

lineStyles = {'-',':','--','-.'};

n = 1000;
sQuantity = 'physicalDose';
if sum(strcmp(fieldnames(result),'RBExDose')) > 0 && ~strcmp(pln.bioOptimization,'none')
    sQuantity = 'RBExDose';
end

dvhPoints = linspace(0,max(result.(sQuantity)(:))*1.05,n);
dvh       = NaN * ones(1,n);

for i = 1:numOfVois
    if cst{i,5}.Visible
        indices     = cst{i,4}{1};
        numOfVoxels = numel(indices);
        doseInVoi   = result.(sQuantity)(indices);   

        % fprintf('%3d %20s - Mean dose = %5.2f Gy +/- %5.2f Gy (Max dose = %5.2f Gy, Min dose = %5.2f Gy)\n', ...
        %     cst{i,1},cst{i,2},mean(doseInVoi),std(doseInVoi),max(doseInVoi),min(doseInVoi))

        for j = 1:n
            dvh(j) = sum(doseInVoi > dvhPoints(j));
        end

        dvh = dvh ./ numOfVoxels * 100;

        subplot(211),plot(dvhPoints,dvh,'LineWidth',4,'Color',colorMx(i,:), ...
            'LineStyle',lineStyles{lineStyleIndicator},'DisplayName',cst{i,2});hold on
    end
end

fontSizeValue = 14;
myLegend = legend('show','location','NorthEast');
set(myLegend,'FontSize',10,'Interpreter','none');
legend boxoff


ylim([0 110]);
xlim([0 1.2*max(dvhPoints)]);
set(gca,'YTick',0:20:120)

grid on,grid minor
box(gca,'on');
set(gca,'LineWidth',1.5,'FontSize',fontSizeValue);
ylabel('Volume [%]','FontSize',fontSizeValue)

if strcmp(sQuantity,'physicalDose');
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
res = matRad_calcQualityIndicators(result,cst,pln);

set(table,'ColumnName',fieldnames(res.QI));
set(table,'Data',(squeeze(struct2cell(res.QI)))');

