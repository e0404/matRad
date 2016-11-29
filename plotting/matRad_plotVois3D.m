function patches = matRad_plotVois3D(axesHandle,ct,cst)
%MATRAD_PLOTVOIS3D Summary of this function goes here
%   Detailed explanation goes here

%Calculate the meshgrid
xCoord = ct.resolution.x * (1:ct.cubeDim(1));
yCoord = ct.resolution.y * (1:ct.cubeDim(2));
zCoord = ct.resolution.z * (1:ct.cubeDim(3));

[xMesh,yMesh,zMesh] = meshgrid(xCoord,yCoord,zCoord);

axes(axesHandle);
wasHold = ishold();

hold(axesHandle,'on');

numVois = size(cst,1);

voiColors = colorcube(numVois);
patches = cell(numVois,1);

for voiIx = 1:numVois
    if ~strcmp(cst{voiIx,3},'IGNORED')
        voiCube = zeros(ct.cubeDim);
        voiCube(cst{voiIx,4}{1}) = 1;
        patches{voiIx} = matRad_plotVoi3D(axesHandle,xCoord,yCoord,zCoord,voiCube,voiColors(voiIx,:),0.3);
    end
end

if ~wasHold
    hold(axesHandle,'off');
end


end

