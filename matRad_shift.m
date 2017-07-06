function [projCoord,idx,targetPoint, sourcePoint] = matRad_shift(initIx,dim,sourcePoint_bev,targetPoint_bev,isoCenter, res, Dx, Dz, rotMat)
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
targetPoint_bev(1) = targetPoint_bev(1) + Dx;
targetPoint_bev(3) = targetPoint_bev(3) + Dz;
sourcePoint_bev(1) = sourcePoint_bev(1) + Dz;
sourcePoint_bev(3) = sourcePoint_bev(3) + Dz;

% rotate in world coord sys and shift by isocenter
targetPoint = targetPoint_bev * rotMat' + isoCenter;
sourcePoint = sourcePoint_bev * rotMat' + isoCenter;

% convert index in coordinates
[coord(:,2),coord(:,1),coord(:,3)] = ind2sub(dim,initIx);
coord = bsxfun(@times,coord,res);

% distance of points Bvec from the projection of points coord onto the line between A 
A2Bnorm = (targetPoint-sourcePoint)/norm(targetPoint-sourcePoint);
pointToEndDist = bsxfun(@minus,coord,targetPoint)*A2Bnorm';

% add translation to the extreme of the ray, according to spherical coord,
% in order to obtain the coord of the projected points
% signvec = sign(Bvec-Avec);
projCoord = bsxfun(@plus,targetPoint,pointToEndDist*A2Bnorm);

% round to voxel coords
D = round(bsxfun(@rdivide,projCoord,res));

% delete every point which goes out of the matrix
D( D(:,1)<1 | D(:,1)>dim(1) | D(:,2)<1 | D(:,2)>dim(2) | D(:,3)<1 | D(:,3)>dim(3), :) = [];

% index the found coordinates
 idx = sub2ind(dim,D(:,2),D(:,1),D(:,3));



