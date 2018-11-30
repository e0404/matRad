function matRad_plotDVHBand(nominalDVH, structureStat, doseLabel)
% matRad_plotDVHBand to plot dose volume bands
% 
% call
%   matRad_plotDVHBand(nominalDVH, structureStat, doseLabel)
%
% input
%   nominalDVH:       x axis values
%   structureStat:    lower bound (start of shadowing)
%   doseLabel:        upper bound (end of shadowing)
%
% output
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('doseLabel', 'var') || isempty(doseLabel)
    doseLabel = 'a.u.';
end

figure;
numOfConf = floor(size(structureStat.dvhStat.percDVH,1) / 2);
% DVH
doseGrid = structureStat.dvhStat.mean.doseGrid;

% plot nominal plan
[y, argmin] = cutAtArgmin(nominalDVH.volumePoints);
x = nominalDVH.doseGrid(1:argmin); 
h(1) = plot(x,y,'LineWidth',2, 'Color', 'k', 'DisplayName', 'nominal');
hold on;
% plot mean
[y, argmin] = cutAtArgmin(structureStat.dvhStat.mean.volumePoints);
x = structureStat.dvhStat.mean.doseGrid(1:argmin); 
h(2) = plot(x,y,'--','LineWidth',2, 'Color', 'k', 'DisplayName', '\mu');
% plot dvh confidence bands
% colors
colors = jet(numOfConf);
alphaTrans = 1;

hIx = numel(h);
for j = 1:numOfConf
    hIx = hIx + 1;
    lIx = j;
    hIx = size(structureStat.dvhStat.percDVH,1) - (j-1);
    lowerLimit = structureStat.dvhStat.percDVH(lIx,:);
    upperLimit = structureStat.dvhStat.percDVH(hIx,:);
    confIn = structureStat.percentiles(hIx) - structureStat.percentiles(lIx);
    confName = ['C', num2str(round(confIn * 100,0))];
    h(hIx) = matRad_shadowPlot(doseGrid, lowerLimit, upperLimit, colors(j,:), confName, alphaTrans);
end

ylim([0 100]);
xlabel(doseLabel);

ylabel('Volume [%]');
lh = legend('show','Location','northeastoutside');
uistack(h(2), 'top')
uistack(h(1), 'top')
labels = get(legend(), 'String');
neworder = numel(labels):-1:1;
plots = flipud(get(gca, 'children'));

% Now re-create the legend
legend(plots(neworder), labels(neworder))

drawnow;
hold off;


function [y, argmin] = cutAtArgmin(x)
  [~,argmin] = min(x);
  y          = x(1:argmin);
end

end

