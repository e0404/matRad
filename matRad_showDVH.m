function matRad_showDVH(dvh,cst,pln,lineStyleIndicator)
% matRad dvh visualizaion
% 
% call
%   matRad_showDVH(dvh,cst)
%   matRad_showDVH(dvh,cst,pln)
%   matRad_showDVH(dvh,cst,lineStyleIndicator)
%   matRad_showDVH(dvh,cst,pln,lineStyleIndicator)
%
% input
%   dvh:                result struct from fluence optimization/sequencing
%   cst:                matRad cst struct
%   pln:                (now optional) matRad pln struct,
%                       standard uses Dose [Gy]
%   lineStyleIndicator: (optional) integer (1,2,3,4) to indicate the current linestyle
%                       (hint: use different lineStyles to overlay
%                       different dvhs)
%
% output
%   graphical display of DVH   
%
% References
%   -
%
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

if ~exist('lineStyleIndicator','var') || isempty(lineStyleIndicator)
    lineStyleIndicator = 1;
end

% create new figure and set default line style indicator if not explictly
% specified
hold on;

%reduce cst
visibleIx = cellfun(@(c) c.Visible == 1,cst(:,5));
cstNames = cst(visibleIx,2);
cstInfo = cst(visibleIx,5);
dvh = dvh(visibleIx);

numOfVois = numel(cstNames);
        
%% print the dvh

%try to get colors from cst
try
    colorMx = cellfun(@(c) c.visibleColor,cstInfo,'UniformOutput',false);
    colorMx = cell2mat(colorMx);
catch
    colorMx    = colorcube;
    colorMx    = colorMx(1:floor(64/numOfVois):64,:);
end

lineStyles = {'-',':','--','-.'};

maxDVHvol  = 0;
maxDVHdose = 0;

for i = 1:numOfVois
    % cut off at the first zero value where there is no more signal
    % behind
    ix      = max([1 find(dvh(i).volumePoints>0,1,'last')]);
    currDvh = [dvh(i).doseGrid(1:ix);dvh(i).volumePoints(1:ix)];
    
    plot(currDvh(1,:),currDvh(2,:),'LineWidth',4,'Color',colorMx(i,:), ...
        'LineStyle',lineStyles{lineStyleIndicator},'DisplayName',cstNames{i})
    
    maxDVHvol  = max(maxDVHvol,max(currDvh(2,:)));
    maxDVHdose = max(maxDVHdose,max(currDvh(1,:)));
end

fontSizeValue = 14;
myLegend = legend('show','location','NorthEast');
set(myLegend,'FontSize',10,'Interpreter','none');
legend boxoff

ylim([0 1.1*maxDVHvol]);
xlim([0 1.2*maxDVHdose]);

grid on,grid minor
box(gca,'on');
set(gca,'LineWidth',1.5,'FontSize',fontSizeValue);
ylabel('Volume [%]','FontSize',fontSizeValue)

if exist('pln','var') && ~isempty(pln)
    if strcmp(pln.propOpt.bioOptimization,'none')
        xlabel('Dose [Gy]','FontSize',fontSizeValue);
    else
        xlabel('RBE x Dose [Gy(RBE)]','FontSize',fontSizeValue);
    end
else
     xlabel('Dose [Gy]','FontSize',fontSizeValue);
end
