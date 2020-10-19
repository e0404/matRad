function  matRad_plotAxisLabels(axesHandle,ct,plane,slice,defaultFontSize,tickdist)
% matRad function to plot x and y labels denoting the ct dimensions 
% according to the selected plane
%
% call
%   matRad_plotAxisLabels(axesHandle,ct,plane,slice,defaultFontSize)
%   matRad_plotAxisLabels(axesHandle,ct,plane,slice,defaultFontSize,
%   tickdist)
%
% input
%   axesHandle         handle to axes the slice should be displayed in
%   ct                 matRad ct structure
%   plane              plane view (coronal=1,sagittal=2,axial=3)
%   slice              slice in the selected plane of the 3D cube
%   defaultFontSize    default font size as double value
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
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('defaultFontSize','var') || isempty(defaultFontSize)
    defaultFontSize = 12;
end

if ~exist('tickdist','var') || isempty(tickdist)
    tickdist = 50;
end
%% Set axis labels and plot iso center
if  plane == 3% Axial plane
    if ~isempty(ct.resolution.x) && ~isempty(ct.resolution.y)
        set(axesHandle,'XTick',0:tickdist/ct.resolution.x:1000);
        set(axesHandle,'YTick',0:tickdist/ct.resolution.y:1000);
        set(axesHandle,'XTickLabel',0:tickdist:1000*ct.resolution.x);
        set(axesHandle,'YTickLabel',0:tickdist:1000*ct.resolution.y);   
        xlabel(axesHandle,'x [mm]','FontSize',defaultFontSize)
        ylabel(axesHandle,'y [mm]','FontSize',defaultFontSize)
        title(axesHandle,['axial plane z = ' num2str(ct.resolution.z*slice) ' [mm]'],'FontSize',defaultFontSize)
    else
        xlabel(axesHandle,'x [voxels]','FontSize',defaultFontSize)
        ylabel(axesHandle,'y [voxels]','FontSize',defaultFontSize)
        title(axesHandle,'axial plane','FontSize',defaultFontSize)
    end
elseif plane == 2 % Sagittal plane
    if ~isempty(ct.resolution.y) && ~isempty(ct.resolution.z)
        set(axesHandle,'XTick',0:tickdist/ct.resolution.z:1000)
        set(axesHandle,'YTick',0:tickdist/ct.resolution.y:1000)
        set(axesHandle,'XTickLabel',0:tickdist:1000*ct.resolution.z)
        set(axesHandle,'YTickLabel',0:tickdist:1000*ct.resolution.y)
        xlabel(axesHandle,'z [mm]','FontSize',defaultFontSize);
        ylabel(axesHandle,'y [mm]','FontSize',defaultFontSize);
        title(axesHandle,['sagittal plane x = ' num2str(ct.resolution.y*slice) ' [mm]'],'FontSize',defaultFontSize)
    else
        xlabel(axesHandle,'z [voxels]','FontSize',defaultFontSize)
        ylabel(axesHandle,'y [voxels]','FontSize',defaultFontSize)
        title(axesHandle,'sagittal plane','FontSize',defaultFontSize);
    end
elseif plane == 1 % Coronal plane
    if ~isempty(ct.resolution.x) && ~isempty(ct.resolution.z)
        set(axesHandle,'XTick',0:tickdist/ct.resolution.z:1000)
        set(axesHandle,'YTick',0:tickdist/ct.resolution.x:1000)
        set(axesHandle,'XTickLabel',0:tickdist:1000*ct.resolution.z)
        set(axesHandle,'YTickLabel',0:tickdist:1000*ct.resolution.x)
        xlabel(axesHandle,'z [mm]','FontSize',defaultFontSize)
        ylabel(axesHandle,'x [mm]','FontSize',defaultFontSize)
        title(axesHandle,['coronal plane y = ' num2str(ct.resolution.x*slice) ' [mm]'],'FontSize',defaultFontSize)
    else
        xlabel(axesHandle,'z [voxels]','FontSize',defaultFontSize)
        ylabel(axesHandle,'x [voxels]','FontSize',defaultFontSize)
        title(axesHandle,'coronal plane','FontSize',defaultFontSize)
    end
end

end

