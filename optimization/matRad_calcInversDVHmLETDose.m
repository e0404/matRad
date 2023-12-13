function mLETDose = matRad_calcInversDVHmLETDose(volume,mLETDoseVec)
% matRad inverse DVH mLETDose (mLETDose Volume Histogram) calculation
% 
% call
%   dose = matRad_calcInversDVHmLETDose(volume,mLETDoseVec)
%
% input
%   volume:     rel volume of structure
%   mLETDoseVec:    mLETDose vector of specific structure
%
% output
%   mLETDose:       mLETDose that corresponds to rel volume
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

% sort mLETDose values
mLETDosePoints = sort(mLETDoseVec,'descend');

ix = max([1 ceil(volume*numel(mLETDosePoints))]);

mLETDose = mLETDosePoints(ix);

end
