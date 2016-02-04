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
% Copyright 2016, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the Eclipse Public License 1.0 (EPL-1.0).
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.
%
% You should have received a copy of the EPL-1.0 in the file license.txt
% along with matRad. If not, see <http://opensource.org/licenses/EPL-1.0>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vAlpha = NaN*ones(numel(vRadDepths),1);
vBeta  = NaN*ones(numel(vRadDepths),1);

numOfTissueClass = size(baseData(1).alpha,2);

% range shift
depths = baseData.depths + baseData.offset;

for i = 1:numOfTissueClass
    vAlpha(mTissueClass==i) = interp1(depths,baseData.alpha(:,i),vRadDepths(mTissueClass==i),'linear');
    vBeta(mTissueClass==i)  = interp1(depths,baseData.beta(:,i), vRadDepths(mTissueClass==i),'linear');
end

