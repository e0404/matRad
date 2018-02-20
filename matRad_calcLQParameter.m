function [vAlpha, vBeta] = matRad_calcLQParameter(vRadDepths,mTissueClass,baseData)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad inverse planning wrapper function
% 
% call
%   [vAlpha, vBeta] = matRad_calcLQParameter(vRadDepths,mTissueClass,baseData)
%
% input
%   vRadDepths:     radiological depths of voxels
%   mTissueClass:   tissue classes of voxels
%   baseData:       biological base data
%
% output
%   vAlpha:         alpha values for voxels interpolated from base data
%   vBeta:          beta values for voxels interpolated from base data
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

vAlpha = NaN*ones(numel(vRadDepths),1);
vBeta  = NaN*ones(numel(vRadDepths),1);

numOfTissueClass = size(baseData(1).alpha,2);

% range shift
depths = baseData.depths + baseData.offset;

for i = 1:numOfTissueClass
    mask = mTissueClass == i;
    if any(mask)
        vAlpha(mask) = matRad_interp1(depths,baseData.alpha(:,i),vRadDepths(mask));
        vBeta(mask)  = matRad_interp1(depths,baseData.beta(:,i), vRadDepths(mask));
    end
end

