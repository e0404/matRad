function [projCoord,idx,targetPoint, sourcePoint] = matRad_projectOnComponents(initIx,dim,sourcePoint_bev,targetPoint_bev,isoCenter, res, Dx, Dz, rotMat)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function projects a point on a certain ray and returns both the index
% of the projected point in the reference system of the ct cube and its
% coordinates
%
% call
%   [finalWeight, X1, sigma_sub, radius, posx, posy, numOfSub] = 
%                                  matRad_calcWeights(sigma_ray, n, method)
%
% input
%   initIx:             initial indices of the points
%   cubeDim:            dimension of the ct cube (i.e. ct.cubeDim)
%   sourcePoint_bev:    source point of the ray in bev
%   targetPoint_bev:    target point of the ray in bev
%   isoCenter:          isocenter coordinates
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
% This file is not part of the offical matRad release
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add offset to target and source point in bev

targetPoint_bevVec(:,2) = targetPoint_bev(2).*ones([length(Dx) 1]);
targetPoint_bevVec(:,1) = bsxfun(@plus,Dx,targetPoint_bev(1));
targetPoint_bevVec(:,3) = bsxfun(@plus,Dz',targetPoint_bev(3));
sourcePoint_bevVec(:,2) = sourcePoint_bev(2).*ones([length(Dz) 1]);
sourcePoint_bevVec(:,1) = bsxfun(@plus,Dx,sourcePoint_bev(1));
sourcePoint_bevVec(:,3) = bsxfun(@plus,Dz,sourcePoint_bev(3));

% rotate in world coord sys and shift by isocenter
targetPoint = bsxfun(@plus,targetPoint_bevVec * rotMat', isoCenter);
sourcePoint = bsxfun(@plus,sourcePoint_bevVec * rotMat', isoCenter);

% convert index in coordinates, it does the same job as:
%       [coord(:,2),coord(:,1),coord(:,3)] = ind2sub(dim,initIx);
if max(initIx) > prod(dim) || min(initIx) < 0
    error('Index exceeds matrix dimensions')
else
    v1 = floor(initIx./dim(1)); 
    coord(:,2) = rem(initIx,dim(1)); coord(coord(:,2)==0, 2) = dim(1);
    coord(:,1) = rem(v1,dim(2)) + 1; 
    coord(coord(:,2)==dim(1), 1) = coord(coord(:,2)==dim(1), 1) -1;
    coord(coord(:,1)==0, 1) = dim(2);
    coord(:,3) = floor((v1)./dim(2)) + 1;
end


coord = bsxfun(@times,coord,res);

% distance of points Bvec from the projection of points coord onto the line between A 
A2Bnorm = (targetPoint-sourcePoint)/norm(targetPoint-sourcePoint);

% We just use the first row of targetPoint because the distance from one
% end to the projected point does not change for parallel translations of
% the ray
pointToEndDist = bsxfun(@minus,coord,targetPoint(1,:))*A2Bnorm';

% add translation to the extreme of the ray, according to spherical coord,
% in order to obtain the coord of the projected points
% signvec = sign(Bvec-Avec);
projCoord = bsxfun(@plus,reshape(targetPoint',[1,3,length(Dx)]),pointToEndDist*A2Bnorm);

% round to voxel coords
D = round(bsxfun(@rdivide,projCoord,res));

% delete every point which goes out of the matrix
% D( D(:,1)<1 | D(:,1)>dim(1) | D(:,2)<1 | D(:,2)>dim(2) | D(:,3)<1 | D(:,3)>dim(3), :) = [];

% index the found coordinates, it does the same thing as:
%       idx = sub2ind(dim,D(:,2),D(:,1),D(:,3));
 idx = D(:,2) + (D(:,1)-1)*dim(1) + (D(:,3)-1)*dim(1)*dim(2);



