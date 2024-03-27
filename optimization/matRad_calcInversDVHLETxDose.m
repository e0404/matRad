function LETxDose = matRad_calcInversDVHLETxDose(volume,LETxDoseVec)
% matRad inverse DVH LETxDose (LETxDose Volume Histogram) calculation
% 
% call
%   dose = matRad_calcInversDVHLETxDose(volume,LETxDoseVec)
%
% input
%   volume:         rel volume of structure
%   LETxDoseVec:    LETxDose vector of specific structure
%
% output
%   LETxDose:       LETxDose that corresponds to rel volume
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

% sort LETxDose values
LETxDosePoints = sort(LETxDoseVec,'descend');

ix = max([1 ceil(volume*numel(LETxDosePoints))]);

LETxDose = LETxDosePoints(ix);

end
