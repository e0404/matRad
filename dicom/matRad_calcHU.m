function ct = matRad_calcHU(ct)
% matRad function to calculate Hounsfield units from a dicom ct 
% that originally uses intensity values
%
% call
%   ct = matRad_calcHU(ct)
%
% input
%   ct: unprocessed dicom ct data which are stored as intensity values (IV)
%
%                      HU = IV * slope + intercept
%
% output
%   ct:                 ct struct with cube with HU
%
% References
%   -
%
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

for i = 1:ct.numOfCtScen
    ct.cubeHU{i} = double(ct.cubeIV{i}) * double(ct.dicomInfo.RescaleSlope) + double(ct.dicomInfo.RescaleIntercept);
end

ct = rmfield(ct,'cubeIV');

end
