function [dvfHandle] = matRad_plotDvfSlice(axesHandle,dvfCube,cubeIdx,plane,slice,colorLine)
% matRad function that generates the plot for the CT in the GUI
% The function can also be used in personal matlab figures by passing the
% corresponding axes handle
%
% call
%   [ctHandle,cMap,window] = matRad_plotCtSlice(axesHandle,ctCube,cubeIdx,plane,slice)
%   [ctHandle,cMap,window] = matRad_plotCtSlice(axesHandle,ctCube,cubeIdx,plane,slice,cMap)
%   [ctHandle,cMap,window] = matRad_plotCtSlice(axesHandle,ctCube,cubeIdx,plane,slice,window)
%   [ctHandle,cMap,window] = matRad_plotCtSlice(axesHandle,ctCube,cubeIdx,plane,slice,cMap,window)
%
% input
%   axesHandle  handle to axes the slice should be displayed in
%   ctCube      the cell of ct cubes
%   cubeIdx     Index of the desired cube in the ct struct
%   plane       plane view (coronal=1,sagittal=2,axial=3)
%   slice       slice in the selected plane of the 3D cube
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

%Use default color line
if nargin < 6 || isempty(colorLine)
    colorLine = [1 1 1];
end

[~,xDim,yDim,zDim] = size(dvfCube{1});

%Prepare the slice and convert it to uint8
if plane == 1 % Coronal plane
    [mX,mZ] = meshgrid(1:zDim,1:xDim);
    xVectorField = squeeze(dvfCube{cubeIdx}(1,:,slice,:));
    zVectorField = squeeze(dvfCube{cubeIdx}(3,:,slice,:));
    dvfHandle=quiver(mX,mZ,zVectorField,xVectorField,'color',colorLine,'Parent',axesHandle);
elseif plane == 2 % sagittal plane
    [mY,mZ] = meshgrid(1:zDim,1:yDim);
    yVectorField = squeeze(dvfCube{cubeIdx}(2,:,slice,:));
    zVectorField = squeeze(dvfCube{cubeIdx}(3,:,slice,:));
    dvfHandle=quiver(mY,mZ,zVectorField,yVectorField,'color',colorLine,'Parent',axesHandle);
elseif plane == 3 % Axial plane
    [mX,mY] = meshgrid(1:yDim,1:xDim);
    xVectorField = squeeze(dvfCube{cubeIdx}(1,:,:,slice));
    yVectorField = squeeze(dvfCube{cubeIdx}(2,:,:,slice));
    dvfHandle=quiver(mX,mY,yVectorField,xVectorField,'color',colorLine,'Parent',axesHandle);
else
    matRad_cfg.dispError('Invalid plane ''%d'' selected for visualization!',plane);
end

end

