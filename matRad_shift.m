function [D0,idx,target, source] = matRad_shift(initIx,dim,source,target,iso, res, ang, Dx, Dz, rotMat)
% This function project a point on a certain ray and return the index of
% the projected point in the reference system of the ct cube
%
% The inputs are:
%   - initIx : initial indeces of the points
%   - cubeDim : dimension of the ct cube (i.e. ct.cubeDim)
%   - source : source point of the ray
%   - target : target point of the ray
%   - iso : isocenter coordinates
%   - res : resolution
%   - ang : row vector with gantry 
%   - Dx : displacement on x axis
%   - Dz : displacement on z axis
% 
% The outputs are:
%   - idx : projected indeces
%   - idx_shift : indeces of the ray translated on the secondary sub-beam


% shift for isocenter and new ray
target(1) = target(1) + Dx;
target(3) = target(3) + Dz;
source(1) = source(1) + Dz;
source(3) = source(3) + Dz;

target = target * rotMat + iso;
source = source * rotMat + iso;

% convert index in coordinates
[coord(:,2),coord(:,1),coord(:,3)] = ind2sub(dim,initIx);
coord = coord .* repmat(res,size(coord,1),1);

% this vectors are target and source coord repeated for every row
Avec = repmat(source,length(initIx),1);
Bvec = repmat(target,length(initIx),1);

% % calculates the distance between the point and nearest one on the ray; and
% % the distance between the nearest point on the ray and one extreme of the
% % ray
% perpDist =  sqrt(sum(cross(coord-Avec,coord-Bvec).^2,2))./sqrt(sum((Bvec-Avec).^2,2));
% otherDist = sqrt( sqrt(sum((coord-Bvec).^2,2)).^2 - perpDist.^2 );
% 
% % angles in spherical coordinates
% theta = asin(sqrt(sum((Avec(:,3)-Bvec(:,3)).^2,2))./sqrt(sum((Avec-Bvec).^2,2)));
% phi = acos(sqrt(sum((Avec(:,1)-Bvec(:,1)).^2,2))./sqrt(sum((Avec(:,1:2)-Bvec(:,1:2)).^2,2)));

pointToEndDist = dot((coord-Bvec),(Avec-Bvec),2)./repmat(norm(target-source),[size(coord,1),1]);

% add translation to the extreme of the ray, according to spherical coord,
% in order to obtain the coord of the projected points
% signvec = sign(Bvec-Avec);
D0 = Bvec + (Avec-Bvec)./norm(target-source).*repmat(pointToEndDist,[1,3]);

D = round(D0./repmat(res,[size(D0,1),1]));
% D = D./ bsxfun(@times,ones(size(D)),res);

% delete every point which goes out of the matrix
D( D(:,1)<1 | D(:,1)>dim(1) | D(:,2)<1 | D(:,2)>dim(2) | D(:,3)<1 | D(:,3)>dim(3), :) = [];

% index the found coordinates
idx = sub2ind(dim,D(:,2),D(:,1),D(:,3));

% Shift all indeces. This gives me the indices of the dose points shifted
% on the new ray
% idx_shift = sub2ind(dim, coord(:,2)+(Dx*sind(ang(1))),...
%     coord(:,1)+round(Dx*cosd(ang(1))),...
%     coord(:,3)+round(Dz*cosd(ang(2))) );


