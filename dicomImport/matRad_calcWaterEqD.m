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
%   ct:                 ct struct with cube with relative _electron_ densities 
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


%% conversion from IV to HU
ctHU = double(ct.cube) * double(ct.dicomInfo.RescaleSlope) + double(ct.dicomInfo.RescaleIntercept);

%% conversion from HU to water equivalent density
load hlutDefault.mat; % load LUT

% Manual adjustments if ct data is corrupt. If some values are out of range
% of the LUT, then these values are adjusted.
if max(ctHU(:)) > max(hlut(:,1))
    warning('projecting out of range HU values');
    ctHU(ctHU>max(hlut(:,1))) = max(hlut(:,1));
end
if min(ctHU(:)) < min(hlut(:,1))
    warning('projecting out of range HU values');
    ctHU(ctHU<min(hlut(:,1))) = min(hlut(:,1));
end

% interpolate HU to relative electron density based on lookup table
ct.cube = interp1(hlut(:,1),hlut(:,2),double(ctHU));

% save hlut
ct.hlut = hlut;

end
