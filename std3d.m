load('slabPhantom1.mat')

kernelsize = 5;

imagesc(ct.cube{1}(:,:,80));
% ct.cube{1}(isnan(ct.cube{1})) = 0;
mu5kernel = ones(kernelsize,kernelsize,kernelsize) ./ kernelsize^3;


mu5  = convn(ct.cube{1}, mu5kernel, 'same');
std5 = convn(ct.cube{1}.^2, mu5kernel, 'same') - mu5.^2;

imagesc(std5(:,:,80))
hold on 
% contour(ct.cube{1}(:,:,round(ct.cubeDim(3)/2)),3,'color','white');
% imagesc(mu5(:,:,80));
