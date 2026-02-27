function matRad_showDVH(dvh,cst,varargin)
% matRad dvh visualizaion
% 
% call
%   matRad_showDVH(dvh,cst)
%   matRad_showDVH(dvh,cst,pln)
%   matRad_showDVH(dvh,cst,Name,Value)
%   matRad_showDVH(dvh,cst,pln,Name,Value)
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

%% Parse input
p = inputParser;

p.addRequired('dvh',@isstruct);
p.addRequired('cst',@iscell);
p.addOptional('pln',[],@(x) isstruct(x) || isempty(x));
p.addParameter('axesHandle',[],@isgraphics);
p.addParameter('plotLegend',true,@(x) isscalar(x) && islogical(x));
%p.addParameter('plotObjectives',false,@(x) isscalar(x) && islogical(x));
p.addParameter('LineWidth',2.5,@(x) isscalar(x) && x > 0);
p.CaseSensitive = false;
p.KeepUnmatched = true;

p.parse(dvh,cst,varargin{:});

axesHandle = p.Results.axesHandle;
dvh = p.Results.dvh;
cst = p.Results.cst;
pln = p.Results.pln;

plotLegend = p.Results.plotLegend;
%plotObjectives = p.Results.plotObjectives;

%{
if plotObjectives && isempty(pln)
    matRad_cfg.dispWarning('Plotting objectives requries pln struct! Disabling.');
    plotObjectives = false;
end
%}

lineWidth = p.Results.LineWidth;

%Get unmatched arguments
fields = fieldnames(p.Unmatched);
values = struct2cell(p.Unmatched);
unmatchedPlotArguments = [fields';values'];


%% Preprocessing
if isempty(axesHandle)
    axesHandle = axes(...
        figure('Color',matRad_cfg.gui.backgroundColor),...
        'Color',matRad_cfg.gui.elementColor,...
        'XColor',matRad_cfg.gui.textColor,...
        'YColor',matRad_cfg.gui.textColor,...
        'GridColor',matRad_cfg.gui.textColor,...
        'MinorGridColor',matRad_cfg.gui.backgroundColor);
end

%Hold for plotting multiple DVHs
hold(axesHandle,'on');

hFig = get(axesHandle,'Parent');
h = findobj(hFig,'type','legend');

if ~isempty(h) 
    if ~plotLegend    
        h.AutoUpdate = false;
    else
        h.AutoUpdate = true;
    end
end

%reduce cst
visibleIx = cellfun(@(c) c.Visible == 1,cst(:,5));
cstNames = cst(visibleIx,2);
cstInfo = cst(visibleIx,5);
%cstObjectives = cst(visibleIx,6);
dvh = dvh(visibleIx);

numOfVois = numel(cstNames);
      
%try to get colors from cst
try
    colorMx = cellfun(@(c) c.visibleColor,cstInfo,'UniformOutput',false);
    colorMx = cell2mat(colorMx);
catch
    colorMx    = colorcube;
    colorMx    = colorMx(1:floor(64/numOfVois):64,:);
end

maxDVHvol  = 0;
maxDVHdose = 0;

%% print the dvhs
for i = 1:numOfVois
    % cut off at the first zero value where there is no more signal
    % behind
    ix      = max([1 find(dvh(i).volumePoints>0,1,'last')]);
    currDvh = [dvh(i).doseGrid(1:ix);dvh(i).volumePoints(1:ix)];
    
    if plotLegend
        plot(axesHandle,currDvh(1,:),currDvh(2,:),'LineWidth',lineWidth,'Color',colorMx(i,:),'DisplayName',cstNames{i},unmatchedPlotArguments{:});
    else
        plot(axesHandle,currDvh(1,:),currDvh(2,:),'LineWidth',lineWidth,'Color',colorMx(i,:),unmatchedPlotArguments{:});
    end
    
    maxDVHvol  = max(maxDVHvol,max(currDvh(2,:)));
    maxDVHdose = max(maxDVHdose,max(currDvh(1,:)));
end


%% Legend and limits
fontSizeValue = matRad_cfg.gui.fontSize;
if plotLegend
    myLegend = legend(axesHandle,'location','NorthEast','FontSize',matRad_cfg.gui.fontSize,'TextColor',matRad_cfg.gui.textColor,'Interpreter','none');
    legend(axesHandle,'boxoff');
    legend(axesHandle,'show');
end

ylim(axesHandle,[0 1.1*maxDVHvol]);
xlim(axesHandle,[0 1.2*maxDVHdose]);

grid(axesHandle,'on'),grid(axesHandle,'minor')
box(axesHandle,'on'); %box(gca,'on');
set(axesHandle,'LineWidth',1.5,'FontSize',fontSizeValue); %set(gca,'LineWidth',1.5,'FontSize',fontSizeValue);
ylabel(axesHandle,'Volume [%]','FontSize',fontSizeValue)

if ~isempty(pln)
    if ~isfield(pln,'bioModel')
        pln.bioModel = 'none';
    end

    pln.bioModel = matRad_BiologicalModel.validate(pln.bioModel,pln.radiationMode);
    
    if strcmp(pln.bioModel.model,'none')
        xlabel(axesHandle,'Dose [Gy]','FontSize',fontSizeValue);
    else
        xlabel(axesHandle,'RBE x Dose [Gy(RBE)]','FontSize',fontSizeValue);
    end
else
    xlabel('Dose [Gy]','FontSize',fontSizeValue);
end
hold(axesHandle,'off');
