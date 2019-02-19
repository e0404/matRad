function [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,thresh,alpha,contourColorMap,...
                                                                          doseColorMap,doseWindow,doseIsoLevels,voiSelection,colorBarLabel,boolPlotLegend,varargin)

% matRad tool function to directly plot a complete slice of a ct with dose
% including contours and isolines.
%
% call
%   [cMapHandle,hDose,hCt,hContour] = matRad_plotSliceWrapper(axesHandle,ct,cst,dose,plane,slice,thresh,alpha,contourColorMap,doseColorMap,doseWindow,doseIsoLevels)
%
% input (required)
%   axesHandle      handle to axes the slice should be displayed in
%   ct              matRad ct struct
%   cubeIdx         Index of the desired cube in the ct struct
%   dose            dose cube
%   plane           plane view (coronal=1,sagittal=2,axial=3)
%   slice           slice in the selected plane of the 3D cube
%
% input (optional / empty)
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
%   varargin        additional input parameters that are passed on to
%                   individual plotting functions (e.g. 'LineWidth',1.5)
%   
%
% output
%   hCMap       handle to the colormap
%   hDose       handle to the dose plot
%   hCt         handle to the ct plot
%   hContour    handle to the contour plot
%   hIsoDose    handle to iso dose contours
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

% Handle the argument list
if ~exist('thresh','var') || isempty(thresh)
    thresh = [];
end
if ~exist('alpha','var') || isempty(alpha)
    alpha = [];
end
if ~exist('contourColorMap','var') || isempty(contourColorMap)
   contourColorMap = [];
end
if ~exist('doseColorMap','var') || isempty(doseColorMap)
   doseColorMap = [];
end
if ~exist('doseWindow','var') || isempty(doseWindow)
   doseWindow = [];
end
if ~exist('doseIsoLevels','var') || isempty(doseIsoLevels)
   doseIsoLevels = [];
end

if ~exist('voiSelection','var') || isempty(voiSelection)
   voiSelection = [];
end

if ~exist('colorBarLabel','var') || isempty(colorBarLabel)
   colorBarLabel = [];
end

if ~exist('boolPlotLegend','var') || isempty(boolPlotLegend)
   boolPlotLegend = false;
end

if ~exist('cst','var') || isempty(cst)
   cst = [];
end

set(axesHandle,'YDir','Reverse');
% plot ct slice
hCt = matRad_plotCtSlice(axesHandle,ct.cubeHU,cubeIdx,plane,slice); 
hold on;

% plot dose
[hDose,doseColorMap,doseWindow] = matRad_plotDoseSlice(axesHandle,dose,plane,slice,thresh,alpha,doseColorMap,doseWindow);

% plot iso dose lines
if ~isempty(doseIsoLevels)
    hIsoDose = matRad_plotIsoDoseLines(axesHandle,dose,[],doseIsoLevels,false,plane,slice,doseColorMap,doseWindow,varargin{:});
    hold on;
else
    hIsoDose = [];
end

%plot VOI contours
if  ~isempty(cst)
    hContour = matRad_plotVoiContourSlice(axesHandle,cst,ct.cubeHU,cubeIdx,voiSelection,plane,slice,contourColorMap,varargin{:});

if boolPlotLegend
   visibleOnSlice = (~cellfun(@isempty,hContour));
   ixLegend = find(voiSelection);
   hContourTmp    = cellfun(@(X) X(1),hContour(visibleOnSlice),'UniformOutput',false);
   if ~isempty(voiSelection)
       hLegend        =  legend(axesHandle,[hContourTmp{:}],[cst(ixLegend(visibleOnSlice),2)],'AutoUpdate','off');
   else
       hLegend        =  legend(axesHandle,[hContourTmp{:}],[cst(visibleOnSlice,2)],'AutoUpdate','off');
   end
   set(hLegend,'Box','Off');
   set(hLegend,'TextColor',[1 1 1]);
   set(hLegend,'FontSize',12);

end
else
    hContour = [];
end

axis(axesHandle,'tight');
set(axesHandle,'xtick',[],'ytick',[]);
colormap(axesHandle,doseColorMap);

matRad_plotAxisLabels(axesHandle,ct,plane,slice,[])

% set axis ratio

ratios = [1/ct.resolution.x 1/ct.resolution.y 1/ct.resolution.z];
   
set(axesHandle,'DataAspectRatioMode','manual');
if plane == 1 
      res = [ratios(3) ratios(2)]./max([ratios(3) ratios(2)]);  
      set(axesHandle,'DataAspectRatio',[res 1])
elseif plane == 2 % sagittal plane
      res = [ratios(3) ratios(1)]./max([ratios(3) ratios(1)]);  
      set(axesHandle,'DataAspectRatio',[res 1]) 
elseif  plane == 3 % Axial plane
      res = [ratios(2) ratios(1)]./max([ratios(2) ratios(1)]);  
      set(axesHandle,'DataAspectRatio',[res 1])
end

hCMap = matRad_plotColorbar(axesHandle,doseColorMap,doseWindow,'Location','EastOutside');
if ~isempty(colorBarLabel)
    set(get(hCMap,'YLabel'),'String', colorBarLabel,'FontSize',14);
end

end

