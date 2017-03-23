cube1=dose_3mm;
cubex=dose_3mm;

figure
imagesc(cube1(:,:,36))

cut1 = any(any(cube1,2),3);
cut2 = any(any(cube1,1),3);
cut3 = any(any(cube1,1),2);

cubex( ~cut1, : , : ) = [];  %rows
cubex( :, ~cut2 , : ) = [];  %columns
cubex( :, :, ~cut3) = [];  %columns

figure
imagesc(cubex(:,:,36))

cubey = zeros(size(cube1));
cubey(cut1,cut2,cut3) = cubex;

figure
imagesc(cubey(:,:,36))
