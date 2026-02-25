function [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,thresh,alpha,contourColorMap,...
                                                                          doseColorMap,doseWindow,doseIsoLevels,voiSelection,colorBarLabel,boolPlotLegend,varargin)
% matRad tool function to directly plot a complete slice of a ct with dose
% including contours and isolines.
%
% call
% [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice)
% [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,thresh)
% [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,alpha)
% [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,contourColorMap)
% [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,doseColorMap)
% [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,doseWindow)
% [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,doseIsoLevels)
%               ...
% [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,thresh,alpha,contourColorMap,...
%                                                                          doseColorMap,doseWindow,doseIsoLevels,voiSelection,colorBarLabel,boolPlotLegend,...)
%
% input (required)
%   axesHandle      handle to axes the slice should be displayed in
%   ct              matRad ct struct
%   cst             matRad cst struct
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

matRad_cfg = MatRad_Config.instance();

matRad_cfg.dispDeprecationWarning('Deprecation warning: matRad_plotSliceWrapper is deprecated. Using matRad_plot_Slice instead');

[hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSlice(ct, 'axesHandle', axesHandle, 'cst', cst, 'cubeIdx', cubeIdx, 'dose', dose, 'plane', plane, 'slice', slice,'thresh', thresh, 'alpha', alpha, 'contourColorMap', contourColorMap, 'doseColorMap', doseColorMap, 'doseWindow', doseWindow, 'doseIsoLevels', doseIsoLevels, 'voiSelection', voiSelection, 'colorBarLabel', colorBarLabel, 'boolPlotLegend', boolPlotLegend, 'others', varargin);

end

