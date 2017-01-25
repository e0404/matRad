function patches = matRad_plotVois3D(axesHandle,ct,cst,selection,cMap)
%MATRAD_PLOTVOIS3D Summary of this function goes here
%   Detailed explanation goes here

if size(cst,2) < 8
    cst = matRad_computeAllVoiSurfaces(ct,cst);
end

%Use default colormap?
if nargin < 5 || isempty(cMap)
    cMap = colorcube(size(cst,1));
end

if nargin < 4 || isempty(selection) || numel(selection) ~= size(cst,1)
    selection = logical(ones(size(cst,1),1));
end

cMapScale = size(cMap,1)-1;

%determine colors
voiColors = cMap(round(linspace(1,cMapScale,size(cst,1))),:);

axes(axesHandle);
wasHold = ishold();

hold(axesHandle,'on');

numVois = size(cst,1);

patches = cell(0);

for voiIx = 1:numVois
    if selection(voiIx) && ~strcmp(cst{voiIx,3},'IGNORED')
        patches{voiIx} = patch(axesHandle,cst{voiIx,8}{1},'VertexNormals',cst{voiIx,8}{2},'FaceColor',voiColors(voiIx,:),'EdgeColor','none','FaceAlpha',0.3);
    end
end

if ~wasHold
    hold(axesHandle,'off');
end


end

