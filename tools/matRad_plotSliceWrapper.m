function [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,thresh,alpha,contourColorMap,doseColorMap,colorMapLabel,doseWindow,doseIsoLevels)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%
% output
%   hCMap       handle to the colormap
%   hDose       handle to the dose plot
%   hCt         handle to the ct plot
%   hContour    handle to the contour plot
%   hIsoDose    handle to iso dose contours
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

%Handle the argument list
if ~exist('tresh','var') || isempty(thresh)
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
if ~exist('colorMapLabel','var') || isempty(colorMapLabel)
   colorMapLabel = '';
end
if ~exist('disableInterpretor','var') || isempty(disableInterpretor)
    disableInterpretor = true;
end


set(axesHandle,'YDir','Reverse');
hCt = matRad_plotCtSlice(axesHandle,ct.cube,cubeIdx,plane,slice); 
hold on;

hContour = matRad_plotVoiContourSlice(axesHandle,cst,ct.cube,cubeIdx,[],plane,slice,contourColorMap);
[hDose,doseColorMap,doseWindow] = matRad_plotDoseSlice(axesHandle,dose,plane,slice,thresh,alpha,doseColorMap,doseWindow);
if ~isempty(doseIsoLevels)
    hIsoDose = matRad_plotIsoDoseLines(axesHandle,dose,[],doseIsoLevels,false,plane,slice,doseColorMap,doseWindow);
else
    hIsoDose = [];
end
axis(axesHandle,'tight');
set(axesHandle,'xtick',[],'ytick',[]);
daspect(axesHandle,[1 1 1]);
colormap(doseColorMap);
hCMap = matRad_plotColorbar(axesHandle,doseColorMap,doseWindow,'Location','EastOutside');
end
