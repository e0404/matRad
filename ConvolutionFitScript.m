load totalSigmaFitDataWedge100.mat
load PHANTOM_wedge_entrance_100mm.mat
% load slab_10_ent_sigmaMap.mat

% figure
% subplot(1,2,1);
% imagesc(mcDose(:,:,25));
% subplot(1,2,2);
% imagesc(anaDose(:,:,25));

res = [ct.resolution.y, ct.resolution.z];

kernelSize = min(size(anaDose,2), size(anaDose,3)) + (mod(min(size(anaDose,2), size(anaDose,3)),2) -  1);
sigmaMap = zeros(size(anaDose));

for i = 1:size(anaDose,1)
    i
    tmpDoseSlice = reshape(anaDose(i,:,:),size(anaDose,2),size(anaDose,3));
    tmpDoseSlice = padarray(tmpDoseSlice,[floor(kernelSize/2), floor(kernelSize/2)],0,'both');
    
    for j = floor(kernelSize/2) + 1:floor(kernelSize/2) + size(anaDose,2)
        for k = floor(kernelSize/2) + 1:floor(kernelSize/2) + size(anaDose,3)
            ixDoseI = j - floor(kernelSize/2);
            ixDoseJ = k - floor(kernelSize/2);
            
%             diff = 1e10;
%             for sigma = 0:0.1:12
%                 tmpKernel   = matRad_create2dimGaussKernel(kernelSize, sigma, res);
%                 tmpConv     = tmpDoseSlice(j - floor(kernelSize/2):j + floor(kernelSize/2), k - floor(kernelSize/2):k + floor(kernelSize/2));
%                 tmp = tmpConv .* tmpKernel;
%                 tmpDiff = abs(sum(tmp(:)) - mcDose(i, ixDoseI, ixDoseJ));
%                 if tmpDiff < diff
%                     diff = tmpDiff;
%                     sigmaMap(i, ixDoseI, ixDoseJ) = sigma;
%                 end
%             end

            f = @(sigma) abs(sum(tmpDoseSlice(j - floor(kernelSize/2):j + floor(kernelSize/2), k - floor(kernelSize/2):k + floor(kernelSize/2)) ...
                            .* matRad_create2dimGaussKernel(kernelSize, sigma, res),'all') - mcDose(i, ixDoseI, ixDoseJ));

            sigmaMap(i, ixDoseI, ixDoseJ) = fminbnd(f, 0, 18);
            
        end
    end
end
anaConvDose = matRad_sliceConvnFilter(ct, anaDose,sigmaMap);
imagesc(sigmaMap(:,:,25))
hold on
contour(mcDose(:,:,round(ct.cubeDim(3)/2)),linspace(0,2e-3,10),'color','black');
colorbar;
pbaspect([ct.cubeDim(2) ct.cubeDim(1) ct.cubeDim(3)]);















