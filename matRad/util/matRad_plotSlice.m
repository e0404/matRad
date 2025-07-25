function [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSlice(ct, varargin)
% matRad tool function to directly plot a complete slice of a ct with dose
% optionally including contours and isolines
%
% call
% [] = matRad_plotSlice(ct, dose, varargin)
%
% input (required)
%   ct              matRad ct struct
%
% input (optional/empty) to be called as Name-value pair arguments:
%   dose            dose cube
%   axesHandle      handle to axes the slice should be displayed in
%   cst             matRad cst struct
%   cubeIdx         Index of the desired cube in the ct struct
%   plane           plane view (coronal=1,sagittal=2,axial=3)
%   slice           slice in the selected plane of the 3D cube
%   thresh          threshold for display of dose values
%   alpha           alpha value for the dose overlay
%   contourColorMap colormap for the VOI contours
%   doseColorMap    colormap for the dose
%   doseWindow      dose value window
%   doseIsoLevels   levels defining the isodose contours
%   voiSelection    logicals defining the current selection of contours
%                   that should be plotted. Can be set to [] to plot
%                   all non-ignored contours.
%   colorBarLabel   string defining the yLabel of the colorBar
%   boolPlotLegend  boolean if legend should be plottet or not
%   showCt          boolean if CT slice should be displayed or not
%   varargin        Additional MATLAB Line or Text Properties (e.g. 'LineWidth', 'FontSize', etc.)
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2025 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultDose             = [];
defaultCst              = [];
defaultSlice            = floor(min(ct.cubeDim)./2);
defaultAxesHandle       = gca;
defaultCubeIdx          = 1;
defaultPlane            = 1;
defaultDoseWindow       = [];
defaultThresh           = [];
defaultAlpha            = [];
defaultDoseColorMap     = jet;
defaultDoseIsoLevels    = [];
defaultVOIselection     = [];
defaultContourColorMap  = [];
defaultBoolPlotLegend   = false;
defaultColorBarLabel    = [];
defaultShowCt           = true;

isDose              = @(x) isnumeric(x) && all(size(x) == ct.cubeDim);
isSlice             = @(x) x>=1 && x<=max(ct.cubeDim) && floor(x)==x;
isAxes              = @(x) strcmp(get(x, 'type'), 'axes');
isCubeIdx           = @(x) isscalar(x);
isPlane             = @(x) isscalar(x) && (sum(x==[1, 2, 3])==1);
isDoseWindow        = @(x) (length(x) == 2 && isvector(x));
isThresh            = @(x) (isscalar(x) && (x>=0) && (x<=1)) || isempty(x);
isAlpha             = @(x) isscalar(x) && (x>=0) && (x<=1) || isempty(x);
isDoseColorMap      = @(x) (isnumeric(x) && (size(x, 2)==3) &&  all(x(:) >= 0) && all(x(:) <= 1)) || isempty(x);
isDoseIsoLevels     = @(x) isnumeric(x) && isvector(x)|| isempty(x);
isVOIselection      = @(x) isnumeric(x) || isempty(x); %all(x(:)==1 | x(:)==0) || isempty(x);
isContourColorMap   = @(x) (isnumeric(x) && (size(x, 2)==3) && size(x, 1)>=2 && all(x(:) >= 0) && all(x(:) <= 1)) || isempty(x);
isBoolPlotLegend    = @(x) x==0 || x ==1;
isColorBarLabel     = @(x) isstring(x) || ischar(x) || isempty(x);
isShowCt            = @(x) isscalar(x) && (x==0) || (x==1);

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'ct')

addParameter(p, 'dose', defaultDose, isDose)
addParameter(p, 'cst', defaultCst)
addParameter(p, 'slice', defaultSlice, isSlice)
addParameter(p, 'axesHandle', defaultAxesHandle, isAxes)
addParameter(p, 'cubeIdx', defaultCubeIdx, isCubeIdx)
addParameter(p, 'plane', defaultPlane, isPlane)
addParameter(p, 'doseWindow', defaultDoseWindow, isDoseWindow)
addParameter(p, 'thresh', defaultThresh, isThresh)
addParameter(p, 'alpha', defaultAlpha, isAlpha)
addParameter(p, 'doseColorMap', defaultDoseColorMap, isDoseColorMap)
addParameter(p, 'doseIsoLevels', defaultDoseIsoLevels, isDoseIsoLevels)
addParameter(p, 'voiSelection', defaultVOIselection, isVOIselection)
addParameter(p, 'contourColorMap', defaultContourColorMap, isContourColorMap)
addParameter(p, 'boolPlotLegend', defaultBoolPlotLegend, isBoolPlotLegend)
addParameter(p, 'colorBarLabel', defaultColorBarLabel, isColorBarLabel)
addParameter(p, 'showCt', defaultShowCt, isShowCt)

parse(p, ct, varargin{:});

%% Unmatched properties
% General properties
lineFieldNames  = fieldnames(set(line));
textFieldNames  = fieldnames(set(text));
axesFieldNames  = fieldnames(set(axes));
% Filter line properties from Unmatched
unmParamNames   = fieldnames(p.Unmatched);
lineFields      = unmParamNames(ismember(unmParamNames, lineFieldNames));
lineValues      = struct2cell(p.Unmatched);
lineValues      = lineValues(ismember(unmParamNames, lineFieldNames));
lineVarargin    = reshape([lineFields, lineValues]', 1, []);
% Filter text properties from Unmatched
textFields      = unmParamNames(ismember(unmParamNames, textFieldNames));
textValues      = struct2cell(p.Unmatched);
textValues      = textValues(ismember(unmParamNames, textFieldNames));
textVarargin    = reshape([textFields, textValues]', 1, []);
%
axesFields      = unmParamNames(ismember(unmParamNames, axesFieldNames));
axesValues      = struct2cell(p.Unmatched);
axesValues      = axesValues(ismember(unmParamNames, axesFieldNames));
axesVarargin    = reshape([axesFields, axesValues]', 1, []);


%% Plot ct slice
matRad_cfg = MatRad_Config.instance();

% Flip axes direction
set(gca,'XTick',[],'YTick',[]);
set(gca,'YDir','Reverse');
% plot ct slice
if p.Results.showCt
    hCt = matRad_plotCtSlice(p.Results.axesHandle,p.Results.ct.cubeHU,p.Results.cubeIdx,p.Results.plane,p.Results.slice, [], []);
end
axis(gca, 'off');

hold on;

%% Plot dose
if ~isempty(p.Results.dose)
    doseWindow = [min(p.Results.dose(:)) max(p.Results.dose(:))];
    if ~isempty(p.Results.doseWindow) && p.Results.doseWindow(2) - p.Results.doseWindow(1) <= 0
        doseWindow = p.Results.doseWindow;
    end        
    [hDose,doseColorMap,doseWindow] = matRad_plotDoseSlice(p.Results.axesHandle, p.Results.dose, p.Results.plane, p.Results.slice, p.Results.thresh, p.Results.alpha, p.Results.doseColorMap, doseWindow);
    hold on;
    
    %% Plot iso dose lines
    hIsoDose = [];
    if ~isempty(p.Results.doseIsoLevels)
        hIsoDose = matRad_plotIsoDoseLines(p.Results.axesHandle,p.Results.dose,[],p.Results.doseIsoLevels,false,p.Results.plane,p.Results.slice,p.Results.doseColorMap,p.Results.doseWindow, lineVarargin{:});
    end

    %% Set Colorbar
    hCMap = matRad_plotColorbar(p.Results.axesHandle,doseColorMap,doseWindow,'Location','EastOutside');
    set(hCMap,'Color',matRad_cfg.gui.textColor);
    if ~isempty(p.Results.colorBarLabel)
        set(get(hCMap,'YLabel'),'String', p.Results.colorBarLabel,'FontSize',matRad_cfg.gui.fontSize);
    end
    set(get(hCMap,'YLabel'),'String', p.Results.colorBarLabel, textVarargin{:});
end
%% Plot VOI contours & Legend

if  ~isempty(p.Results.cst)
    [hContour,~] = matRad_plotVoiContourSlice(p.Results.axesHandle, p.Results.cst, p.Results.ct, p.Results.cubeIdx, p.Results.voiSelection, p.Results.plane, p.Results.slice, p.Results.contourColorMap, lineVarargin{:});

    if p.Results.boolPlotLegend
        visibleOnSlice = (~cellfun(@isempty,hContour));        
        hContourTmp    = cellfun(@(X) X(1),hContour(visibleOnSlice),'UniformOutput',false);
        voiSelection = visibleOnSlice;
        if ~isempty(p.Results.voiSelection)
            voiSelection = visibleOnSlice(find(p.Results.voiSelection));
        end
        hLegend        =  legend(p.Results.axesHandle,[hContourTmp{:}],[p.Results.cst(voiSelection,2)],'AutoUpdate','off','TextColor',matRad_cfg.gui.textColor);
        set(hLegend,'Box','On');
        set(hLegend,'TextColor',matRad_cfg.gui.textColor);
        if ~isempty(textVarargin)
            set(hLegend, textVarargin{:});
        else
            set(hLegend,'FontSize',matRad_cfg.gui.fontSize);
        end
    end
else
    hContour = [];
end

%% Adjust axes
axis(p.Results.axesHandle,'tight');
set(p.Results.axesHandle,'xtick',[],'ytick',[]);
%colormap(p.Results.axesHandle,p.Results.doseColorMap);
fontSize = [];
if isfield(p.Unmatched, 'FontSize')
    fontSize = p.Unmatched.FontSize;
end
matRad_plotAxisLabels(p.Results.axesHandle,p.Results.ct,p.Results.plane,p.Results.slice, fontSize, [])

% Set axis ratio.
ratios = [1/p.Results.ct.resolution.x 1/p.Results.ct.resolution.y 1/p.Results.ct.resolution.z];

set(p.Results.axesHandle,'DataAspectRatioMode','manual');
if p.Results.plane == 1
    res = [ratios(3) ratios(2)]./max([ratios(3) ratios(2)]);
    set(p.Results.axesHandle,'DataAspectRatio',[res 1])
elseif p.Results.plane == 2 % sagittal plane
    res = [ratios(3) ratios(1)]./max([ratios(3) ratios(1)]);
    set(p.Results.axesHandle,'DataAspectRatio',[res 1])
elseif  p.Results.plane == 3 % Axial plane
    res = [ratios(2) ratios(1)]./max([ratios(2) ratios(1)]);
    set(p.Results.axesHandle,'DataAspectRatio',[res 1])
end

%% Set text properties
if ~isempty(textVarargin)
    set(p.Results.axesHandle, textVarargin{:})
    set(p.Results.axesHandle.Title, textVarargin{:})
end

if ~isempty(axesVarargin)
    index = cellfun(@(x) strcmp(x,'Title'), axesVarargin, 'UniformOutput', 1);
    index = find(index);
    p.Results.axesHandle.Title.String = axesVarargin{index+1};
end

if ~exist('hCMap', 'var')
    hCMap = [];
end
if ~exist('hDose', 'var')
    hDose = [];
end
if ~exist('hCt', 'var')
    hCt = [];
end
if ~exist('hContour', 'var')
    hContour = [];
end
if ~exist('hIsoDose', 'var')
    hIsoDose = [];
end

end