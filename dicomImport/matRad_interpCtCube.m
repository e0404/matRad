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
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
