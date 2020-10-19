function patches = matRad_plotVois3D(axesHandle,ct,cst,selection,cMap)
% matRad function that plots 3D structures of the volumes of interest
% If the 3D-data is not stored in the CT, it will be commputed on the fly.
%
% call
%   patches = matRad_plotVois3D(axesHandle,ct,cst,selection)
%   patches = matRad_plotVois3D(axesHandle,ct,cst,selection,cMap)
%
% input
%   axesHandle  handle to axes the structures should be displayed in
%   ct          matRad ct struct which contains resolution
%   cst         matRad cst struct
%   selection   logicals defining the current selection of contours
%               that should be plotted. Can be set to [] to plot
%               all non-ignored contours.
%   cMap        optional argument defining the colormap, default are the
%               colors stored in the cst
%
% output
%   patches     patch objects created by the matlab 3D visualization
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

if size(cst,2) < 8
    cst = matRad_computeAllVoiSurfaces(ct,cst);
end

%Use stored colors or colormap?
if nargin < 5 || isempty(cMap)
    cMapScale = size(cMap,1)-1;
    %determine colors
    voiColors = cMap(round(linspace(1,cMapScale,size(cst,1))),:);
else
    for i = 1:size(cst,1)
      voiColors(i,:) = cst{i,5}.visibleColor;
    end
end

if nargin < 4 || isempty(selection) || numel(selection) ~= size(cst,1)
    selection = logical(ones(size(cst,1),1));
end

cMapScale = size(cMap,1)-1;

axes(axesHandle);
wasHold = ishold();

hold(axesHandle,'on');

numVois = size(cst,1);

patches = cell(0);

for voiIx = 1:numVois
    if selection(voiIx) && ~strcmp(cst{voiIx,3},'IGNORED')
        patches{voiIx} = patch(cst{voiIx,8}{1},'VertexNormals',cst{voiIx,8}{2},'FaceColor',voiColors(voiIx,:),'EdgeColor','none','FaceAlpha',0.4,'Parent',axesHandle);
    end
end

if ~wasHold
    hold(axesHandle,'off');
end


end
