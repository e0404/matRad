function dose = matRad_sliceConvnFilter(ct, standardDose, cStd)

kernelSize = min(ct.cubeDim(2), ct.cubeDim(3)) + (mod(min(ct.cubeDim(2), ct.cubeDim(3)),2) -  1);
res = [ct.resolution.y, ct.resolution.z];
% kernel = matRad_create2dimGaussKernel(kernelSize, sigma, res);

newDose = zeros(size(standardDose));
for ixSlice = 1:size(standardDose,1)
    doseSlice = reshape(standardDose(ixSlice,:,:),ct.cubeDim(2),ct.cubeDim(3));
    doseSlice = padarray(doseSlice,[floor(kernelSize/2), floor(kernelSize/2)],0,'both');
    
    for ii = floor(kernelSize/2) + 1:floor(kernelSize/2) + size(newDose,2)
        for ij = floor(kernelSize/2) + 1:floor(kernelSize/2) + size(newDose,3)
            ixDoseI = ii - floor(kernelSize/2);
            ixDoseJ = ij - floor(kernelSize/2);
            sigma = cStd(ixSlice, ixDoseI, ixDoseJ);
            if (sigma > 0.001)
                tmpKernel   = matRad_create2dimGaussKernel(kernelSize, sigma, res);
                tmpConv     = doseSlice(ii - floor(kernelSize/2):ii + floor(kernelSize/2), ij - floor(kernelSize/2):ij + floor(kernelSize/2));
                tmp = tmpConv .* tmpKernel;
                newDose(ixSlice,ixDoseI, ixDoseJ) = sum(tmp(:));
            else
                newDose(ixSlice,ixDoseI, ixDoseJ) = standardDose(ixSlice,ixDoseI, ixDoseJ);
            end
        end
    end
end

dose = newDose;