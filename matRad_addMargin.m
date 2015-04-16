function [ mVOIEnlarged ] = matRad_addMargin(mVOI,vResolution,sMargin,DiaElem)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad steering information generation
% 
% call
%   stf = matRad_generateStf(ct,cst,pln,visMode)
%
% input
%   mVOI:           image stack in dimensions of X x Y x Z holding ones for
%                   object and zeros otherwise 
%   sMargin:        margin in mm 
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
    DiaElem = false;
elseif nargin <3
    error('not enough input parameters specified for matRad_addMargin');
end

% get number of voxels which should be added in each dimension
VoxelMargins =round(sMargin./vResolution);
mVOIEnlarged=mVOI;
NewIdx = [];

    for Cnt = 1:max(VoxelMargins)
        
        % for multiple loops just consider just added margin
        NewIdx = setdiff(find(mVOIEnlarged),NewIdx);
        [xCoord, yCoord, zCoord] = ind2sub(size(mVOIEnlarged),NewIdx);

        for i=1:numel(xCoord)

            if VoxelMargins(1)>=Cnt
                dx=1;
            else
                dx=0;
            end

            if VoxelMargins(2)>=Cnt
                dy=1;
            else
                dy=0;
            end

            if VoxelMargins(3)>=Cnt
                dz=1;
            else
                dz=0;
            end

            for j= -1:1:1
                    mVOIEnlarged(xCoord(i)   ,yCoord(i)   ,zCoord(i)+dz*j)=1;
                    
                    mVOIEnlarged(xCoord(i)+dx,yCoord(i)   ,zCoord(i)+dz*j)=1;
                    mVOIEnlarged(xCoord(i)-dx,yCoord(i)   ,zCoord(i)+dz*j)=1;
                    mVOIEnlarged(xCoord(i)   ,yCoord(i)+dy,zCoord(i)+dz*j)=1;
                    mVOIEnlarged(xCoord(i)   ,yCoord(i)-dy,zCoord(i)+dz*j)=1;
                if DiaElem
                    mVOIEnlarged(xCoord(i)+dx,yCoord(i)+dy,zCoord(i)+dz*j)=1;
                    mVOIEnlarged(xCoord(i)+dx,yCoord(i)-dy,zCoord(i)+dz*j)=1;
                    mVOIEnlarged(xCoord(i)-dx,yCoord(i)+dy,zCoord(i)+dz*j)=1;
                    mVOIEnlarged(xCoord(i)-dx,yCoord(i)-dy,zCoord(i)+dz*j)=1;
                end
            end
            

        end
    end


end

