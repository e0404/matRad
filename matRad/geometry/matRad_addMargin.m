function mVOIEnlarged = matRad_addMargin(mVOI,cst,vResolution,vMargin,bDiaElem)
% matRad add margin function
% 
% call
%   mVOIEnlarged = matRad_addMargin(mVOI,cst,vResolution,vMargin,bDiaElem)
%
% input
%   mVOI:           image stack in dimensions of X x Y x Z holding ones for
%                   object and zeros otherwise 
%   cst:            matRad cst struct
%   vResolution     ct resolution
%   vMargin:        margin in mm 
%   bDiaElem        if true 26-connectivity is used otherwise 6-connectivity
%
% output
%   mVOIEnlarged:   enlarged VOI
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

if nargin == 4
    bDiaElem = false;
elseif nargin < 4
    error('not enough input parameters specified for matRad_addMargin');
end

% generate voi cube for patient surface/patient skin
idx             = [cst{:,4}];
voiSurfaceIdx   = unique(vertcat(idx{:}));

% get number of voxels which should be added in each dimension
voxelMargins = round([vMargin.x vMargin.y vMargin.z]./[vResolution.x vResolution.y vResolution.z]);
mVOIEnlarged = mVOI;
newIdx = [];

[yUpperLim,xUpperLim,zUpperLim] = size(mVOI);

for cnt = 1:max(voxelMargins)

    % for multiple loops consider just added margin
    newIdx = setdiff(find(mVOIEnlarged),newIdx);
    [yCoord, xCoord, zCoord] = ind2sub(size(mVOIEnlarged),newIdx);

    % find indices on border and take out
    borderIx = xCoord==1 | xCoord==xUpperLim | ...
               yCoord==1 | yCoord==yUpperLim | ...
               zCoord==1 | zCoord==zUpperLim;

    xCoord(borderIx) = [];
    yCoord(borderIx) = [];
    zCoord(borderIx) = [];

    dx = voxelMargins(1)>=cnt;
    dy = voxelMargins(2)>=cnt;
    dz = voxelMargins(3)>=cnt;

    for i = -1:1
        for j = -1:1
            for k = -1:1

                if (abs(i)+abs(j)+abs(k) == 0) || (~bDiaElem && abs(i)+abs(j)+abs(k) > 1) % skip if diagonal elements not wanted or zero offset
                    continue;
                end

                newIx = (yCoord+i*dy) + (xCoord+j*dx-1)*size(mVOIEnlarged,1) + ...
                        (zCoord+k*dz-1)*size(mVOIEnlarged,1)*size(mVOIEnlarged,2);
                    
                % check if new indices are part of voiSurfaceIdx
                bWithinPatient = ismember(newIx,voiSurfaceIdx);
                
                mVOIEnlarged(newIx(bWithinPatient)) = 1;

            end
        end
    end


end

