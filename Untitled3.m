% s = load('mri');
% mriVolume = squeeze(s.D);
% sizeIn = size(mriVolume);
% hFigOriginal = figure;
% hAxOriginal  = axes;
% slice(double(mriVolume),sizeIn(2)/2,sizeIn(1)/2,sizeIn(3)/2);
% grid on, shading interp, colormap gray

cubo = ct.cube{1};
dim = size(cubo);

theta = stf.gantryAngle;
phi = stf.couchAngle;

t1 = [ 1 0 0 0 
    0 cosd(phi) sind(phi) 0
    0 -sind(phi) cosd(phi) 0
    0 0 0 1];

tform = affine3d(t1);

% rocubo = imwarp(cubo,tform);
% 
% if size(rocubo,1) ~= dim(1)
% else if size(rocubo,1) == dim(1) + 1
%         
%     c = round((size(rocubo,1)-dim(1))/2);
%     rocubo(1:c,:,:) = [];
%     rocubo(end-c+1:end,:,:) = [];
    
rocubo = imrotate(cubo,40,'crop');

rocubo2 = permute(rocubo,[1 3 2]);

rocubo2 = imrotate(rocubo2,-10,'crop');

rocubo = permute(rocubo2,[1 3 2]);



t2 = [ cosd(theta) -sind(theta) 0 0
    sind(theta) cosd(theta) 0 0
    0 0 1 0
    0 0 0 1];

tform = affine3d(t2);

rocubo = imwarp(rocubo,tform);