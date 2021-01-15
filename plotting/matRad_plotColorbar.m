function cBarHandle = matRad_plotColorbar(axesHandle,cMap,window,varargin)
% matRad function wrapper for plotting the colorbar 
% This is necessary since the rgb colors are manually mapped within the ct 
% and the dose plotting, and MATLAB attaches colorbars to axes.
%
% call
%   cBarHandle = matRad_plotColorbar(axesHandle,cMap,window,varargin)
%
% input
%   axesHandle  handle to axes the colorbar will be attached to
%   cMap        corresponding colormap
%   window      colormap window (corresponds to clim)
%   varargin    additional key-value pairs that will be forwarded to the
%               MATLAB colorbar(__) call
%
% output
%   cBarHandle  handle of the colorbar object
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

v=version;

colormap(axesHandle,cMap);
caxis(window);

if str2double(v(1:3))>=8.5
    cBarHandle = colorbar(axesHandle,varargin{:});
else
    cBarHandle = colorbar('peer',axesHandle,varargin{:});
end





end

