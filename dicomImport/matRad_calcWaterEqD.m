function ct = matRad_calcWaterEqD(ct)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate the equivalent densities from a 
% dicom ct that originally uses intensity values
%
% call
%   ct = matRad_calcWaterEqD(ct, slope, intercept)
%
% input
%   ct:                 unprocessed dicom ct data which are stored as
%                       intensity values (IV)
%   solpe:              parameter for linear conversion into HU
%   intercept:          parameter for linear conversion into HU
%
%                      HU = IV * slope + intercept
%
% base data
%   HU2waterEqT.mat:    look up table to convert from HU to relative 
%                       electron densities. Note that this is just example
%                       data. If you want to make precise computaitons for
%                       your scanner you need to replace this lookup table
%
% output
%   ctEqD:              ct cube with relative _electron_ densities 
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% using a look up table (LUT) for HU to waterEqThickness, a new cube is 
% generated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% conversion from IV to HU

% reshaping the cube produces real values, but for the conversion using the
% LUT integer values are needed
ctHU = int32(ct.cube * ct.dicomInfo.RescaleSlope + ct.dicomInfo.RescaleIntercept);

%% conversion from HU to water equivalent density
load hlutDefault.mat; % load LUT

% Manual adjustments if ct data is corrupt. If some values are out of range
% of the LUT, then these values are adjusted.
if sum(ctHU > max(hlut(:,1)) | ctHU < min(hlut(:,1)))>0
    warning('projecting out of range HU values');
    ctHU(ctHU<min(hlut(:,1))) = min(hlut(:,1));
    ctHU(ctHU>max(hlut(:,1))) = max(hlut(:,1));
end

% interpolate HU to relative electron density based on lookup table
ct.cube = interp1(hlut(:,1),hlut(:,2),double(ctHU));

% save hlut
ct.hlut = hlut;

end