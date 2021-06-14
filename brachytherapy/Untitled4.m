[x,y] = meshgrid(-10:0.1:10,-10:0.1:10);  % just to get some x and y values
z = sqrt(x.^2 + y.^2);
surf(x,y,z)
shading interp
set(gca,'cameraposition',[0 0 180])  % this essentially turns a 3d surface plot into a 2d plot
colormap(flipud(jet))