function interpCt = matRad_interpCtCube(origCt, origCtInfo, resolution)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call interpCube = interp3dCube(original ct-Cube, origRes, newRes) to 
% change the dimensions of the ct-cube. origRes and newRes are row-vectors
% with the distance of the voxels in mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coordsOfFirstPixel = [origCtInfo.ImagePositionPatient];

% set up grid vectors
x = coordsOfFirstPixel(1,1) + origCtInfo(1).PixelSpacing(1)*double([0:origCtInfo(1).Rows-1]);
y = coordsOfFirstPixel(2,1) + origCtInfo(1).PixelSpacing(2)*double([0:origCtInfo(1).Columns-1]');
z = coordsOfFirstPixel(3,:);

xi = coordsOfFirstPixel(1,1):resolution(1):(coordsOfFirstPixel(1,1)+origCtInfo(1).PixelSpacing(1)*double(origCtInfo(1).Rows-1));
yi = [coordsOfFirstPixel(2,1):resolution(2):(coordsOfFirstPixel(2,1)+origCtInfo(1).PixelSpacing(2)*double(origCtInfo(1).Columns-1))]';
zi = coordsOfFirstPixel(3,1):resolution(3): coordsOfFirstPixel(3,end);

% interpolate
interpCt.cube = interp3(x,y,z,origCt,xi,yi,zi);

% some meta information
interpCt.resolution = resolution;

interpCt.x = xi;
interpCt.y = yi';
interpCt.z = zi;
