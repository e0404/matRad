function hpatch = matRad_plotIsoDose3D(axesHandle,xMesh,yMesh,zMesh,doseCube,isoLevels,cMap,window,alpha)
% matRad function that plots isolines in 3d
%
% call
%   hpatch = matRad_plotIsoDose3D(axesHandle,xMesh,yMesh,zMesh,doseCube,isoLevels,cMap,window,alpha)
%
% input
%   axesHandle  handle to axes the slice should be displayed in
%   x/y/zMesh   meshs
%   doseCube    dose cube
%   isoLevels   levels for computation
%   cMap        colormap
%   window      window for dose display
%   alpha       transparency
%
% output
%   hpatch: handle to the patch object
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

if nargin < 9
    alpha = 0.1;
end

if nargin < 8
    window = [min(doseCube(:)) max(doseCube(:))];
end

if nargin < 7
    cMap = jet(64);
end

if nargin < 6
    isoLevels = linspace(window(1),window(2),5);
end

cMapScale = size(cMap,1) - 1;

isoColorLevel = uint8(cMapScale*(isoLevels - window(1))./(window(2)-window(1)));

%This circumenvents a bug in Octave when the index in the image hase the maximum value of uint8
if matRad_cfg.isOctave
	isoColorLevel(isoColorLevel == 255) = 254;
end

colors = squeeze(ind2rgb(isoColorLevel,cMap));

for isoValIx = 1:numel(isoLevels)
    hpatch = patch(axesHandle,isosurface(xMesh,yMesh,zMesh,doseCube,isoLevels(isoValIx)));
    reducepatch(hpatch);
    isonormals(xMesh,yMesh,zMesh,doseCube,hpatch);
    hpatch.FaceColor = colors(isoValIx,:);
    hpatch.EdgeColor = 'none';
    hpatch.FaceAlpha = alpha;
end
end

