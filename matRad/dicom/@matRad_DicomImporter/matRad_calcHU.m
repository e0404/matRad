function obj = matRad_calcHU(obj)
% matRad function to calculate Hounsfield units from a dicom ct 
% that originally uses intensity values
%
% In your object, there must be a property that contains unprocessed 
% dicom ct data which are stored as intensity values (IV)
%
% Output - ct structure with cube with HU
%
% HU = IV * slope + intercept
%
% call
%   obj = matRad_calcHU(obj)
%
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:obj.ct.numOfCtScen
    obj.ct.cubeHU{i} = double(obj.ct.cubeIV{i}) * double(obj.ct.dicomInfo.RescaleSlope) + double(obj.ct.dicomInfo.RescaleIntercept);
end

obj.ct = rmfield(obj.ct,'cubeIV');

end
