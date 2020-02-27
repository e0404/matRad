clear
load slab_10_ent_sigmaMap.mat
load PHANTOM_hetSlab_entrance_10mm.mat
load radDepthsSlab10.mat

res = [ct.resolution.y, ct.resolution.z];

kernelSize = min(size(radDepths,2), size(radDepths,3)) + (mod(min(size(radDepths,2), size(radDepths,3)),2) -  1);

for i = 1:size(radDepths,1)
    i
    tmpDepthSlice = reshape(radDepths(i,:,:),size(radDepths,2),size(radDepths,3));
    tmpDepthSlice = padarray(tmpDepthSlice,[floor(kernelSize/2), floor(kernelSize/2)],0,'both');
    
    for j = floor(kernelSize/2) + 1:floor(kernelSize/2) + size(radDepths,2)
        for k = floor(kernelSize/2) + 1:floor(kernelSize/2) + size(radDepths,3)
            ixDoseI = j - floor(kernelSize/2);
            ixDoseJ = k - floor(kernelSize/2);
            

            f = @(weight) abs(sqrt( ...
                    sum(tmpDepthSlice(j - floor(kernelSize/2):j + floor(kernelSize/2), k - floor(kernelSize/2):k + floor(kernelSize/2)).^2 ...
                            .* matRad_create2dimGaussKernel(kernelSize, weight, res),'all') - ...
                    sum(tmpDepthSlice(j - floor(kernelSize/2):j + floor(kernelSize/2), k - floor(kernelSize/2):k + floor(kernelSize/2)) ...
                            .* matRad_create2dimGaussKernel(kernelSize, weight, res),'all')^2) - sigmaMap(i, ixDoseI, ixDoseJ));


            weightMap(i, ixDoseI, ixDoseJ) = fminbnd(f, 0, 100);
            
        end
    end
end















