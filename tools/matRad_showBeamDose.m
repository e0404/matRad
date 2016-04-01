function matRad_showBeamDose(dij,w,slice)

tmpDose = reshape(dij.physicalDose{1}*w,dij.dimensions); 
    
subplot(2,2,1)
imagesc(tmpDose(:,:,slice))
axis equal tight
axis([30 140 40 120])
xlabel('voxel index')
ylabel('voxel index')
colorbar
title('sum dose')    

for i = 1:dij.numOfBeams
      
    tmpW = zeros(dij.totalNumOfBixels,1);
    tmpW(dij.beamNum==i) = w(dij.beamNum==i);
    
    tmpDose = reshape(dij.physicalDose{1}*tmpW,dij.dimensions); 
    
    subplot(2,2,i+1)
    imagesc(tmpDose(:,:,slice))
    axis equal tight
    axis([30 140 40 120])
    xlabel('voxel index')
    ylabel('voxel index')
    colorbar
    title(['beam # ' num2str(i)])    
    
end

