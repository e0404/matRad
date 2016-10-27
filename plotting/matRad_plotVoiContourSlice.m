function voiContourHandles = matRad_plotVoiContourSlice(axesHandle,cst,selection,plane,slice,cMap)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function that plots the contours of the segmentations given in cst
%
% call
%   voiContourHandles = matRad_plotVoiContourSlice(axesHandle,cst,plane,slice,cMap)
%
% input
%   axesHandle          handle to axes the slice should be displayed in
%   cst                 cst structure
%   selection           logicals defining the current selection of contours
%                       that should be plotted. Can be set to [] to plot
%                       all non-ignored contours.
%   plane               plane view (coronal=1,sagittal=2,axial=3)
%   slice               slice in the selected plane of the 3D cube
%   cMap                optional argument defining the colormap, default is
%                       colorcube
%
% output
%   voiContourHandles:  handles of the plotted contours
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

%Use default colormap?
if nargin < 6 || isempty(cMap)
    cMap = colorcube(size(cst,1));
end

if isempty(selection) || numel(selection) ~= size(cst,1)
    selection = logical(ones(size(cst,1),1));
end

cMapScale = size(cMap,1)-1;

%determine colors
colors = cMap(round(linspace(1,cMapScale,size(cst,1))),:);

axes(axesHandle)

voiContourHandles = gobjects(0);

for s = 1:size(cst,1)
    if ~strcmp(cst{s,3},'IGNORED') && selection(s)
        if size(cst,2) >= 7 && ~isempty(cst{s,7})
            % plot precalculated contourc data
            if any(cst{s,7}{slice,plane}(:))
                lower = 1; % lower marks the beginning of a section
                while lower-1 ~= size(cst{s,7}{slice,plane},2);
                    hold on
                    steps = cst{s,7}{slice,plane}(2,lower); % number of elements of current line section
                    voiContourHandles(end+1) = line(cst{s,7}{slice,plane}(1,lower+1:lower+steps),...
                        cst{s,7}{slice,plane}(2,lower+1:lower+steps),...
                        'Color',colors(s,:),'LineWidth',2.0,'Parent',axesHandle);
                    
                    lower = lower+steps+1;
                end
            end
        else
            %TODO: use standard contour function
        end
    end
end

end
