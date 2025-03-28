function  matRad_plotAxisLabels(axesHandle,ct,plane,slice,defaultFontSize,titleString)
% matRad function to plot x and y labels denoting the ct dimensions 
% according to the selected plane
%
% call
%   matRad_plotAxisLabels(axesHandle,ct,plane,slice,defaultFontSize)
%   matRad_plotAxisLabels(axesHandle,ct,plane,slice,defaultFontSize,titleString)
%
% input
%   axesHandle         handle to axes the slice should be displayed in
%   ct                 matRad ct structure
%   plane              plane view (coronal=1,sagittal=2,axial=3)
%   slice              slice in the selected plane of the 3D cube
%   defaultFontSize    default font size as double value
%   titleString        optional custom title to display above plane information
%
% output
%   -
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
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

if ~exist('defaultFontSize','var') || isempty(defaultFontSize)
    defaultFontSize = matRad_cfg.gui.fontSize;
end

if ~exist('titleString','var') || isempty(titleString)
    titleString = '';
end

ct = matRad_getWorldAxes(ct);

%% Set axis labels and plot iso center
if plane == 3 % Axial plane
    if ~isempty(ct.resolution.x) && ~isempty(ct.resolution.y)
        set(axesHandle,'XTick',linspace(0,ct.x(end)-ct.x(1),10)./ct.resolution.x); 
        set(axesHandle,'YTick',linspace(0,ct.y(end)-ct.y(1),10)./ct.resolution.y);
        set(axesHandle,'XTickLabel',round(linspace(ct.x(1),ct.x(end),10)));
        set(axesHandle,'YTickLabel',round(linspace(ct.y(1),ct.y(end),10)));   
        xlabel(axesHandle,'x [mm]','FontSize',defaultFontSize)
        ylabel(axesHandle,'y [mm]','FontSize',defaultFontSize)
        vcoord = matRad_cubeIndex2worldCoords([1,1,slice],ct);
        planeInfo = ['axial plane z = ' num2str(vcoord(3)) ' [mm]'];
    else
        xlabel(axesHandle,'x [voxels]','FontSize',defaultFontSize)
        ylabel(axesHandle,'y [voxels]','FontSize',defaultFontSize)
        planeInfo = 'axial plane';
    end
elseif plane == 2 % Sagittal plane
    if ~isempty(ct.resolution.y) && ~isempty(ct.resolution.z)
        set(axesHandle,'XTick',linspace(0,ct.z(end)-ct.z(1),10)./ct.resolution.z);
        set(axesHandle,'YTick',linspace(0,ct.y(end)-ct.y(1),10)./ct.resolution.y);
        set(axesHandle,'XTickLabel',round(linspace(ct.z(1),ct.z(end),10)));
        set(axesHandle,'YTickLabel',round(linspace(ct.y(1),ct.y(end),10)));
        xlabel(axesHandle,'z [mm]','FontSize',defaultFontSize);
        ylabel(axesHandle,'y [mm]','FontSize',defaultFontSize);
        vcoord = matRad_cubeIndex2worldCoords([1,slice,1],ct);
        planeInfo = ['sagittal plane x = ' num2str(vcoord(1)) ' [mm]'];
    else
        xlabel(axesHandle,'z [voxels]','FontSize',defaultFontSize)
        ylabel(axesHandle,'y [voxels]','FontSize',defaultFontSize)
        planeInfo = 'sagittal plane';
    end
elseif plane == 1 % Coronal plane
    if ~isempty(ct.resolution.x) && ~isempty(ct.resolution.z)
        set(axesHandle,'XTick',linspace(0,ct.z(end)-ct.z(1),10)./ct.resolution.z);
        set(axesHandle,'YTick',linspace(0,ct.x(end)-ct.x(1),10)./ct.resolution.x);
        set(axesHandle,'XTickLabel',round(linspace(ct.z(1),ct.z(end),10)));
        set(axesHandle,'YTickLabel',round(linspace(ct.x(1),ct.x(end),10)));
        xlabel(axesHandle,'z [mm]','FontSize',defaultFontSize)
        ylabel(axesHandle,'x [mm]','FontSize',defaultFontSize)
        vcoord = matRad_cubeIndex2worldCoords([slice,1,1],ct);
        planeInfo = ['coronal plane y = ' num2str(vcoord(2)) ' [mm]'];
    else
        xlabel(axesHandle,'z [voxels]','FontSize',defaultFontSize)
        ylabel(axesHandle,'x [voxels]','FontSize',defaultFontSize)
        planeInfo = 'coronal plane';
    end
end

% Create the combined title if a custom title is provided
if ~isempty(titleString)
    titleText = sprintf('%s\n%s', titleString, planeInfo);
else
    titleText = planeInfo;
end

% Set the title with the combined text
title(axesHandle, titleText, 'FontSize', defaultFontSize, 'Color', matRad_cfg.gui.highlightColor);

%Apply coloring
set(axesHandle,'XColor',matRad_cfg.gui.textColor,'YColor',matRad_cfg.gui.textColor);

end


