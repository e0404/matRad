function hpatch = matRad_plotVoi3D(axesHandle,xMesh,yMesh,zMesh,voiSegmentation,color,alpha)
%MATRAD_PLOTVOI3D Summary of this function goes here
%   Detailed explanation goes here
    
    if nargin < 7
        alpha = 0.5;
    end

    v = smooth3(voiSegmentation);
    hpatch = patch(axesHandle,isosurface(xMesh,yMesh,zMesh,v,0.5));
    reducepatch(hpatch);
    isonormals(xMesh,yMesh,zMesh,v,hpatch);
    hpatch.FaceColor = color;
    hpatch.EdgeColor = 'none';
    hpatch.FaceAlpha = alpha;
end

