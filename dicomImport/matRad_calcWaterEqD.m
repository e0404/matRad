function ctEqD = matRad_calcWaterEqD(ct, slope, intercept)
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
ctHU = int32(ct * slope + intercept);

%% conversion from HU to water equivalent density
load('HU2waterEqD.mat'); % load LUT

% Manual adjustments if ct data is corrupt. If some values are out of range
% of the LUT, then these values are adjusted.
ctHU(ctHU<-1024) = -1024;
ctHU(ctHU>3071) = 3071;

% execute conversion to waterEqT: 
% LUT:| waterEqT  |   HU    |
%     |   0       |  -1024  |
%     |   0       |  -1023  |
%     |  ...      |   ...   |

% the first row corresponds to a HU-value of -1024 -> index = ctHU+1025
ctEqD = reshape(HU2waterEqD(ctHU+1025,1),size(ctHU));

end