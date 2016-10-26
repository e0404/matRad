function ctHandle = matRad_plotCtSlice(axesHandle,ctCube,plane,slice,cMap)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function that generates the plot for the CT in the GUI. The
% function can also be used in personal matlab figures by passing the
% corresponding axes handle
%
% call
%   ctHandle = matRad_plotCtSlice(axesHandle,ctCube,plane,slice,cMap)
%
% input
%   axesHandle  handle to axes the slice should be displayed in
%   ctCube      3D array of the CT to select the slice from
%   plane       plane view (coronal=1,sagittal=2,axial=3)
%   slice       slice in the selected plane of the 3D cube
%   cMap        optional argument defining the colormap, default is bone
%
% output
%   ctHandle:   handle of the plotted CT axes
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
if nargin < 5 || isempty(cMap)
    cMap = bone(64);
end

cMapScale = size(cMap,1) - 1;

maxCT = max(ctCube(:));
if plane == 1 % Coronal plane
    ct_rgb = ind2rgb(uint8(cMapScale*(squeeze(ctCube(slice,:,:)/maxCT))),cMap);
      
elseif plane == 2 % sagittal plane
    ct_rgb = ind2rgb(uint8(cMapScale*(squeeze(ctCube(:,slice,:)/maxCT))),cMap);
       
elseif plane == 3 % Axial plane
    ct_rgb = ind2rgb(uint8(cMapScale*(squeeze(ctCube(:,:,slice)/maxCT))),cMap);
end
axes(axesHandle)
ctHandle = imagesc(ct_rgb);

end

