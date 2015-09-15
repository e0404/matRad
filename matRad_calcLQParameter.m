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

vAlpha = NaN*ones(numel(vRadDepths),1);
vBeta  = NaN*ones(numel(vRadDepths),1);

numOfTissueClass = size(baseData(1).alpha,2);

% range shift
depths = baseData.depths + baseData.offset;
Idx = depths > 0;

for i = 1:numOfTissueClass
    vAlpha(mTissueClass==i) = interp1(depths(Idx),baseData.alpha(Idx,i),vRadDepths(mTissueClass==i),'linear');
    vBeta(mTissueClass==i)  = interp1(depths(Idx),baseData.beta(Idx,i), vRadDepths(mTissueClass==i),'linear');
end

