function interpCt = matRad_interpCtCube(origCt, origCtInfo, resolution)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to interpolate a 3D ct cube to a different resolution
% 
% call
%   interpCt = matRad_interpCtCube(origCt, origCtInfo, resolution)
%
% input
%   origCt:         original CT as matlab 3D array
%   origCtInfo:     meta information about the geometry of the orgiCt cube
%   resolution:     target resolution [mm] in x, y, an z direction for the 
%                   new cube
%
% output
%   interpCt:       interpolated ct cube as matlab 3D array
%
% References
%   -
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


coordsOfFirstPixel = [origCtInfo.ImagePositionPatient];

% set up grid vectors
x = coordsOfFirstPixel(1,1) + origCtInfo(1).PixelSpacing(1)*double([0:origCtInfo(1).Columns-1]);
y = coordsOfFirstPixel(2,1) + origCtInfo(1).PixelSpacing(2)*double([0:origCtInfo(1).Rows-1]');
z = coordsOfFirstPixel(3,:);

xi = coordsOfFirstPixel(1,1):resolution.x:(coordsOfFirstPixel(1,1)+origCtInfo(1).PixelSpacing(1)*double(origCtInfo(1).Rows-1));
yi = [coordsOfFirstPixel(2,1):resolution.y:(coordsOfFirstPixel(2,1)+origCtInfo(1).PixelSpacing(2)*double(origCtInfo(1).Columns-1))]';
zi = coordsOfFirstPixel(3,1):resolution.z: coordsOfFirstPixel(3,end);

% interpolate
interpCt.cube = interp3(x,y,z,origCt,xi,yi,zi);

% some meta information
interpCt.resolution = resolution;

interpCt.x = xi;
interpCt.y = yi';
interpCt.z = zi;
