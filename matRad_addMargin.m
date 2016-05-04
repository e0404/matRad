function mVOIEnlarged = matRad_addMargin(mVOI,cst,vResolution,vMargin,bDiaElem)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
voiSurface      = zeros(size(mVOI));
idx             = [cst{:,4}];
idx             = unique(vertcat(idx{:}));
voiSurface(idx) = 1;
voiSurfaceIdx   = find(voiSurface);

% get number of voxels which should be added in each dimension
voxelMargins = round([vMargin.y vMargin.x vMargin.z]./[vResolution.y vResolution.x vResolution.z]);
mVOIEnlarged = mVOI;
NewIdx = [];

[xUpperLim,yUpperLim,zUpperLim]=size(mVOI);

for Cnt = 1:max(voxelMargins)

    % for multiple loops consider just added margin
    NewIdx = setdiff(find(mVOIEnlarged),NewIdx);
    [xCoord, yCoord, zCoord] = ind2sub(size(mVOIEnlarged),NewIdx);

    % find indices on border and take out
    borderIx = xCoord==1 | xCoord==xUpperLim | ...
               yCoord==1 | yCoord==yUpperLim | ...
               zCoord==1 | zCoord==zUpperLim;

    xCoord(borderIx) = [];
    yCoord(borderIx) = [];
    zCoord(borderIx) = [];

    dx = voxelMargins(1)>=Cnt;
    dy = voxelMargins(2)>=Cnt;
    dz = voxelMargins(3)>=Cnt;

    for i = -1:1
        for j = -1:1
            for k = -1:1

                if (abs(i)+abs(j)+abs(k) == 0) || (~bDiaElem && i+j+k > 1) % skip if diagonal elements not wanted or zero offset
                    continue;
                end

                newIx = (xCoord+i*dx) + (yCoord+j*dy-1)*size(mVOIEnlarged,1) + ...
                        (zCoord+k*dz-1)*size(mVOIEnlarged,1)*size(mVOIEnlarged,2);
                    
                % check if new indices are part of voiSurfaceIdx
                bWithinPatient = ismember(newIx,voiSurfaceIdx);
                
                mVOIEnlarged(newIx(bWithinPatient)) = 1;

            end
        end
    end


end

