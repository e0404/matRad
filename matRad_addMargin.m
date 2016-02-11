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

if nargin == 3
    bDiaElem = false;
elseif nargin < 3
    error('not enough input parameters specified for matRad_addMargin');
end

% generate voi cube for patient surface/patient skin
voiSurface = zeros(size(mVOI));
voiSurface(unique([cell2mat(cst(:,4))])) = 1;
voiSurfaceIdx = find(voiSurface);

% get number of voxels which should be added in each dimension
voxelMargins = round([vMargin.x vMargin.y vMargin.z]./[vResolution.x vResolution.y vResolution.z]);
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

                newIx = sub2ind(size(mVOIEnlarged),xCoord+i*dx,yCoord+j*dy,zCoord+k*dz);
                
                % check if new indices are part of voiSurfaceIdx
                bWithinPatient = ismember(newIx,voiSurfaceIdx);
                
                mVOIEnlarged(newIx(bWithinPatient)) = 1;

            end
        end
    end


end

