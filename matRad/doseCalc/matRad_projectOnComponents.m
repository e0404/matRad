function [projCoord,idx,targetPoint, sourcePoint] = matRad_projectOnComponents(initIx,dim,sourcePoint_bev,targetPoint_bev,isoCenter, res, Dx, Dz, rotMat)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function projects a point on a certain ray and returns both the index
% of the projected point in the reference system of the ct cube and its
% coordinates
%
% call
%   [projCoord,idx,targetPoint, sourcePoint] = 
%           matRad_projectOnComponents(initIx,dim,sourcePoint_bev,targetPoint_bev,isoCenter, res, Dx, Dz, rotMat)
%
% input
%   initIx:             initial indices of the points
%   cubeDim:            dimension of the ct cube (i.e. ct.cubeDim)
%   sourcePoint_bev:    source point of the ray in bev
%   targetPoint_bev:    target point of the ray in bev
%   isoCenter:          isocenter coordinates (in cube system)
%   res:                resolution in x,y, and z direction
%   Dx:                 displacement on x axis
%   Dz:                 displacement on z axis
%
% output
%   idx:        projected indeces
%   projCoord:  projected coordinates
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add offset to target and source point in bev
d = [Dx zeros(size(Dx)) Dz];
targetPoint = (targetPoint_bev + d) * rotMat' + isoCenter; 
sourcePoint = (sourcePoint_bev + d) * rotMat' + isoCenter; 


% convert index in coordinates, it does the same job as:
%       [coord(:,2),coord(:,1),coord(:,3)] = ind2sub(dim,initIx);
%{
if max(initIx) > prod(dim) || min(initIx) < 0
    error('Index exceeds matrix dimensions')
else
    v1 = floor(initIx ./ dim(1)); 
    coord(:,2) = rem(initIx, dim(1)); coord(coord(:,2) == 0, 2) = dim(1);
    coord(:,1) = rem(v1,dim(2)) + 1; 
    coord(coord(:,2) == dim(1), 1) = coord(coord(:,2) == dim(1), 1) - 1;
    coord(coord(:,1) == 0, 1) = dim(2);
    coord(:,3) = floor((v1) ./ dim(2)) + 1;
end
%}

[coord(:,2),coord(:,1),coord(:,3)] = ind2sub(dim,initIx);

%Multiply with resolution
coord = coord.*res;

% distance of points Bvec from the projection of points coord onto the line between A 
A2Bnorm = (targetPoint-sourcePoint)/norm(targetPoint-sourcePoint);

% We just use the first row of targetPoint because the distance from one
% end to the projected point does not change for parallel translations of
% the ray
pointToEndDist = (coord - targetPoint(1,:))*A2Bnorm';

% add translation to the extreme of the ray, according to spherical coord,
% in order to obtain the coord of the projected points
% signvec = sign(Bvec-Avec);
%projCoord = bsxfun(@plus,reshape(targetPoint',[1,3,length(Dx)]),pointToEndDist*A2Bnorm);
projCoord = reshape(targetPoint',[1 3 size(targetPoint,1)]) + pointToEndDist*A2Bnorm;

% round to voxel coords
% D = round(bsxfun(@rdivide,projCoord,res));
D = round(projCoord./res);

% delete every point which goes out of the matrix
% D( D(:,1)<1 | D(:,1)>dim(1) | D(:,2)<1 | D(:,2)>dim(2) | D(:,3)<1 | D(:,3)>dim(3), :) = [];

% index the found coordinates, it does the same thing as:
%       idx = sub2ind(dim,D(:,2),D(:,1),D(:,3));
 idx = D(:,2) + (D(:,1)-1)*dim(1) + (D(:,3)-1)*dim(1)*dim(2);


