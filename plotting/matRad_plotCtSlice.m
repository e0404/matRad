function [ctHandle,cMap,window] = matRad_plotCtSlice(axesHandle,ctCube,cubeIdx,plane,slice,cMap,window)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function that generates the plot for the CT in the GUI. The
% function can also be used in personal matlab figures by passing the
% corresponding axes handle
%
% call
%   [ctHandle,cMap,window] = matRad_plotCtSlice(axesHandle,ct,cubeIdx,plane,slice,cMap,window)
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%Use default colormap?
if nargin < 6 || isempty(cMap)
    cMap = bone(64);
end

if nargin < 7 || isempty(window)
    window = [min(ctCube{cubeIdx}(:)) max(ctCube{cubeIdx}(:))];    
end

cMapScale = size(cMap,1) - 1;

if plane == 1 % Coronal plane
    ct_rgb = ind2rgb(uint8(cMapScale*(squeeze((ctCube{cubeIdx}(slice,:,:)-window(1))/(window(2) - window(1))))),cMap);
      
elseif plane == 2 % sagittal plane
    ct_rgb = ind2rgb(uint8(cMapScale*(squeeze((ctCube{cubeIdx}(:,slice,:)-window(1))/(window(2) - window(1))))),cMap);
       
elseif plane == 3 % Axial plane
    ct_rgb = ind2rgb(uint8(cMapScale*(squeeze((ctCube{cubeIdx}(:,:,slice)-window(1))/(window(2) - window(1))))),cMap);
end
ctHandle = image('CData',ct_rgb,'Parent',axesHandle);

end

