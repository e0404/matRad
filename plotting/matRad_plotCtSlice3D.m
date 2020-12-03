function [ctHandle,cMap,window] = matRad_plotCtSlice3D(axesHandle,ct,cubeIdx,plane,ctSlice,cMap,window)
% matRad function that generates the plot for the CT in the GUI 3D view 
% The function can also be used in personal matlab figures by passing the
% corresponding axes handle
%
% call
%   [ctHandle,cMap,window] = matRad_plotCtSlice3D(axesHandle,ct,cubeIdx,plane,ctSlice,cMap,window)
%
% input
%   axesHandle  handle to 3D axes the slice should be displayed in
%   ct          the ct struct used in matRad
%   cubeIdx     Index of the desired cube in the ct struct
%   plane       plane view (coronal=1,sagittal=2,axial=3)
%   ctSlice     slice in the selected plane of the 3D cube
%   cMap        optional argument defining the colormap, default is bone
%               if you want to use the default map with the window argument
%               you can use an empty array []
%   window      optional argument defining the displayed range. default is
%               [min(ctCube(:)) max(ctCube(:))]
%
% output
%   ctHandle    handle of the plotted CT axes
%   cMap        used colormap (same as input if set)
%   window      used window (same as input if set)
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

%Use default colormap?
if nargin < 6 || isempty(cMap)
    cMap = bone(64);
end

if nargin < 7 || isempty(window)
    window = [min(ct.cubeHU{cubeIdx}(:)) max(ct.cubeHU{cubeIdx}(:))];    
end

cMapScale = size(cMap,1) - 1;

%Create the coordinates
coords{1} = ct.resolution.x * (1:ct.cubeDim(2));
coords{2} = ct.resolution.y * (1:ct.cubeDim(1));
coords{3} = ct.resolution.z * (1:ct.cubeDim(3));
 
% slice plot with surface(...), colormapping can be done by texture
% mapping, this is why we use surface instead of slice
if plane == 1 % Coronal plane
    [xMesh,zMesh] = meshgrid(coords{2},coords{3});
    yMesh = ctSlice*ct.resolution.x*ones(size(xMesh));
    ctIndexed = uint8(cMapScale*(squeeze((ct.cubeHU{cubeIdx}(ctSlice,:,:)-window(1))/(window(2) - window(1)))));
%    ct_rgb = permute(ct_rgb,[2 1 3]);
elseif plane == 2 % sagittal plane
    [yMesh,zMesh] = meshgrid(coords{1},coords{3});
    xMesh = ctSlice*ct.resolution.y*ones(size(yMesh));
    ctIndexed = uint8(cMapScale*(squeeze((ct.cubeHU{cubeIdx}(:,ctSlice,:)-window(1))/(window(2) - window(1)))));
%    ct_rgb = permute(ct_rgb,[2 1 3]);
elseif plane == 3 % Axial plane
    [xMesh,yMesh] = meshgrid(coords{2},coords{1});
    zMesh = ctSlice*ct.resolution.z*ones(size(xMesh));    
    ctIndexed = uint8(cMapScale*(squeeze((ct.cubeHU{cubeIdx}(:,:,ctSlice)-window(1))/(window(2) - window(1)))));
else
	matRad_cfg.dispError('Invalid plane ''%d'' selected for visualization!',plane);
end

%This circumenvents a bug in Octave when the index in the image hase the maximum value of uint8
if matRad_cfg.isOctave
	ctIndexed(ctIndexed == 255) = 254;
end

ct_rgb = ind2rgb(ctIndexed,cMap);

ctHandle = surface('XData',xMesh, 'YData',yMesh, 'ZData',zMesh, ...
        'CData',ct_rgb, 'CDataMapping','direct', ...
        'EdgeColor','none', 'FaceColor','texturemap', 'BackFaceLighting','unlit','FaceLighting','flat','Parent',axesHandle);

%{
% slice plot with slice(...)
% seems technically more reasonable, however manual colormapping not straightforward        
[xMesh,yMesh,zMesh] = meshgrid(coords{:});
slicePlane{1} = [];
slicePlane{2} = [];
slicePlane{3} = [];

slicePlane{plane} = coords{plane}(ctSlice);

ct_rgb = ind2rgb(uint8(cMapScale*((ct.cubeHU{cubeIdx}-window(1))/(window(2) - window(1)))),cMap);

ctHandle = slice(axesHandle,xMesh,yMesh,zMesh,ct.cubeHU{cubeIdx},slicePlane{:});

set(ctHandle,'EdgeColor','none');
%}



end



