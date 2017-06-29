function [TracMat]=matRad_multipleRayTracing(radDepthMat,rayEn,iso,res)

pos = rayEn{2} + repmat(iso,[size(rayEn{2},1) 1]);
pos = pos ./repmat(res,[size(pos,1) 1]);

dim = size(radDepthMat);

pos2 = round(repmat(pos,[1 1 dim(2)]));

x = 1:dim(2);

xx = repmat(x,[size(rayEn{2},1) 1]);

pos2(:,2,:) = xx;

idx = sub2ind(dim, reshape(pos2(:,2,:),1,[]), reshape(pos2(:,1,:),1,[]), reshape(pos2(:,3,:),1,[]));

TracMat = zeros(size(radDepthMat));
TracMat(idx) = radDepthMat(idx);