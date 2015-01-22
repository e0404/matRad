function [ct] = matRad_genCT(result)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function reads a CT from MGH Data set.
% Rescaling of DICOM data to a user-defined intercept, slope, and data type.
% Save the CT with anatomical (LPS) coordinate system .
% Interpolate Head and neck CT.
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

% Load .mat file who has all CT data.
if result == 1
    load('prostate3mmvoxels.mat');
elseif result == 2
    load('liver3mmvoxels.mat');
elseif result == 3
    load('hn3mmvoxels.mat');
elseif result == 4
    load('tg119_core_and_outertarget_3mmvoxels.mat');
end


%% Save CT data on ctCube variable.
ctCube = planC{3}.scanArray(:,:,:);

%% Head and neck interpolation
% Head and neck CT and VOI countours not have same dimensions and scale. To fix
% it, we need to interpolate the CT.

if result == 3

% Save CT in 'V' variable.
V = double(ctCube);

% Save the X and Y dimensions 
[dimX,dimY,~] = size(V);

% Save slices Z position .
z = [planC{3}.scanInfo(:).zValue];

% Create a mesgrid whit same CT dimensions.
[X,Y,Z]=meshgrid(1:dimX,1:dimY,z);

% Create a vector with initial Z position until the last Z position to be
% interpolated.
zi = 64.25:.5:(64.25+66*.5);

% Create a mesgrid whit CT interpolated dimensions.
[XI,YI,ZI]=meshgrid(1:dimX,1:dimY,zi);

% 3D interpolate.
ctCube=interp3(X,Y,Z,V,XI,YI,ZI,'spline');

% The interpolations gives reals numbers, but we need intergers.
ctCube=int64(ctCube);

end

%% Rescaling of DICOM data to a user-defined intercept, slope, and data type.
% CT data are in Intensities values, so it is necessary to transform to
% Hounsfield units.
% For CT data, intercept and slope define the transformation between intensity
% values (IV) and Hounsfield Units (HU): HU = IV * slope + intercept.
% More information: http://goo.gl/KpEX6q

% Define CTAir CTWater.

CTAir = planC{1,3}.uniformScanInfo.CTAir;
CTWater = planC{1,3}.uniformScanInfo.CTWater;

% Rescaling from intensities values to Hounsfield Units. Note: The original CT is saved
% in unit16 format, because intensities values are positive. It is necessary transform to double,
% because in next step (LUT for HU to WE) LUT only accept doubles.

% Here we calculate the intercept and slope from CT data, here is calculate
% the transformation (HU): HU = IV * slope + intercept.
% Below the linear interpolation is obtained using knows HU values of CT for
% Air and Water: CT_Air = -1000 and, CT_Water = 0
ct  = (double(ctCube(:,:,:))-CTAir)*((1000)/(CTWater-CTAir)) - 1000;


%% Convert Hounsfield Units to water equivalent.

% Load HU2waterEqT.mat, contains a look up table (LUT) to transform HU to
% water equivalent units.
load('HU2waterEqT.mat');

% Manual adjustments when ct data corrupt. If some value is out of range for
% LUT table, then this vales are adjusted.
ct(ct<-1024) = -1024;
ct(ct>3071) = 3071;

%Apply LUT conversion, the output of this conversion gives a 1-column
%vector, so it is necessary form again the CT cube. This is done by reshape
%command.
% LUT input are in HU. Note: we need to apply a shift of +1025, because
% LUT started in -1024.
ct = reshape(HU2waterEqT(ct+1025,1),size(ct));

%% Save the CT with anatomical (LPS) coordinate system.
% We use the most important model for coordinate system for medical
% imaging, the original data do not use this system, it is only consistent
% for X and Y axis. It is needed to mirror every slice in Z axis.
% More information about LPS system: http://goo.gl/gPM9YP

% This for cicle, runs over every Z slice reflecting it. For example for a
% Cube who has 90 slices, the slice z=89, is reflected to slice z=2 and
% viceversa.

ct_temp = zeros(size(ct));

for j=1:size(ct,3)
    ct_temp(:,:,size(ct,3)-j+1) = ct(:,:,j);
end

% Saved reflected CT.
ct = ct_temp;