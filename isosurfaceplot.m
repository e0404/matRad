cubo = zeros(dij.dimensions);
% cubo = zeros(size(radDepthCube{1}));
% cubo(V(ix)) = radDepths;
% cubo(V) = radDepthV{1};
% cubo(V) = radDepthCube{1}(V);
% cubo = radDepthCube{1};

cubo(idx) = radDepthsMat(idx);

cubo(isnan(cubo)) = 0;

maxi = max(max(max(cubo)));
isovec = [0.001*maxi, 0.2*maxi, 0.5*maxi, 0.8*maxi, 0.95*maxi];
facevec = [0.1, 0.1, 0.2, 0.4, 1];
colorvec = ['c', 'b', 'g', 'y', 'r'];

for i=1:5
    [f,v] = isosurface(cubo, isovec(i));
    patch('Faces',f,'Vertices',v,'EdgeAlpha',0.0001,'FaceAlpha',facevec(i),'FaceColor',colorvec(i));
end
