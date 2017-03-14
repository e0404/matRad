function hpatch = matRad_plotIsoDose3D(axesHandle,xMesh,yMesh,zMesh,doseCube,isoLevels,cMap,window,alpha)
%MATRAD_PLOTVOI3D Summary of this function goes here
%   Detailed explanation goes here
    
    if nargin < 9
        alpha = 0.1;
    end
    
    if nargin < 8
        window = [min(doseCube(:)) max(doseCube(:))];
    end
    
    if nargin < 7
        cMap = jet(64);
    end
    
    if nargin < 6
        isoLevels = linspace(window(1),window(2),5);
    end
    
    cMapScale = size(cMap,1) - 1;
    isoColorLevel = (isoLevels - window(1))./(window(2)-window(1));
    isoColorLevel(isoColorLevel < 0) = 0;
    isoColorLevel(isoColorLevel > 1) = 0;
    colors = squeeze(ind2rgb(uint8(cMapScale*isoColorLevel),cMap));
    
    
    for isoValIx = 1:numel(isoLevels)
        hpatch = patch(axesHandle,isosurface(xMesh,yMesh,zMesh,doseCube,isoLevels(isoValIx)));
        reducepatch(hpatch);
        isonormals(xMesh,yMesh,zMesh,doseCube,hpatch);
        hpatch.FaceColor = colors(isoValIx,:);
        hpatch.EdgeColor = 'none';
        hpatch.FaceAlpha = alpha;
    end
end

