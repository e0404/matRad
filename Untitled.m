% i have a ray
% convert to coord

tic 
xmin = 80;
xmax = 500;
y0 = 250;
z0 = 21;
rd = 30;
res = [1 1 10];

K = zeros(512,512,39);
dim = size(K);
[x,y,z] = meshgrid(1:size(K,1),1:size(K,2),1:size(K,3));

V = [reshape(x,prod(dim),1), reshape(y,prod(dim),1), reshape(z,prod(dim),1)];

V( V(:,1)<xmin | V(:,1)>xmax | (y0-V(:,2)).^2+(z0*res(3)-V(:,3)*res(3)).^2 > rd^2, :) = [];

idx = sub2ind(dim, V(:,2), V(:,1), V(:,3));
toc