function movementPlotting(ct, slice)

figure

for i = 1:ct.numOfCtScen
    
    im = ct.cube{1}(:,:,slice);
    subplot(121)
    imagesc(im)
    subplot(122)
    imagesc(ct.cube{i}(:,:,slice))
    title(i)
    pause(.2)
    hold on 
    
end

end