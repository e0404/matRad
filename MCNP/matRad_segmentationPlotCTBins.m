function matRad_segmentationPlotCTBins(ct, plane, sliceOfInt)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad visualization of tissue segmentation on CT-data
%
% call
%   matRad_segmentationPlotCTBins(ct, plane, slice)
%
% input
%       ct:         CT-data with HU in ct.cubeHU and linear indices of
%                   segmenation in ct.tissueBin.linIndVol and ct.cubeDim
%       plane:      view options ('axial' (default), 'sagital', 'coronal')
%       slice:      slice you would like to view (default = half of
%                   dimension)
%
% output
%   ...
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 03/2019
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prelude
if ~exist('plane')
    plane = 'axial';
    planeOption = 3;
elseif strcmp(plane, 'sagital')
        planeOption = 2;
elseif strcmp(plane, 'coronal')
        planeOption = 1;
end

if ~exist('sliceOfInt')
    sliceOfInt = round(ct.cubeDim(planeOption)/2);
end


segmentationVolume = zeros(ct.cubeDim);

for matCounter = 1:size({ct.tissueBin.matIndex},2)
segmentationVolume(ct.tissueBin(matCounter).linIndVol) = matCounter -1;
end

% Plot

fig1 = figure;
fig1.Units = 'centimeters';
fig1.Position = [0 0 29.7 21];
ax1 = axes;
if strcmp(plane, 'axial')
    ctPlot = imagesc(squeeze(ct.cubeHU{1}(:,:,sliceOfInt)));
elseif strcmp(plane, 'sagital')
    ctPlot = imagesc(squeeze(ct.cubeHU{1}(:,sliceOfInt,:)));
elseif strcmp(plane, 'coronal')
    ctPlot = imagesc(squeeze(ct.cubeHU{1}(sliceOfInt,:,:)));
end

title('Segmentation of CT-Scan Using Hounsfield Unit Intervals', 'FontSize', 20)


ax2 = axes;
linkaxes([ax1,ax2])
if strcmp(plane, 'axial')
    segmentationPlot = contour(flip(squeeze(segmentationVolume(:,:,sliceOfInt)),1), [1:size({ct.tissueBin.matIndex},2)-1]);
elseif strcmp(plane, 'sagital')
    segmentationPlot = contour(flip(squeeze(segmentationVolume(:,sliceOfInt,:)),1), [1:size({ct.tissueBin.matIndex},2)-1]);
elseif strcmp(plane, 'coronal')
    segmentationPlot = contour(flip(squeeze(segmentationVolume(sliceOfInt,:,:)),1), [1:size({ct.tissueBin.matIndex},2)-1]);
end


ax1.XTick = [];
ax1.YTick = [];
ax1.Units = 'centimeters';

ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
ax2.Units = 'centimeters';


set([ax1,ax2],'Position',[4 1.5 21.2 17]);
cb1 = colorbar(ax1); 
cb1.Units = 'centimeters';
cb1.Position = [3 1.5 .5 17];

cb2 = colorbar(ax2);
cb2.Units = 'centimeters';
cb2.Position = [25.7 1.5 .5 17];
cb2.TickLabels = {ct.tissueBin.name};
cb2.Ticks = [0:size({ct.tissueBin.matIndex},2)-1];

cb1.Label.String = 'Hounsfield Units';

cb1.FontSize = 18;
cb2.FontSize = 18;

colormap(ax1, 'gray')
colormap(ax2, 'autumn')



end