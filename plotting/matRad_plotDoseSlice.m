function [doseHandle,cMap,window] = matRad_plotDoseSlice(axesHandle,doseCube,plane,slice,threshold,alpha,cMap,window)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function that generates a dose plot of a selected slice. The
% function can also be used in personal matlab figures by passing the
% corresponding axes handle
%
% call
%   [doseHandle,cMap,window] = matRad_plotDoseSlice(axesHandle,ctCube,plane,slice,threshold,alpha,cMap,window)
%
% input
%   axesHandle  handle to axes the slice should be displayed in
%   doseCube    3D array of the dose to select the slice from
%   plane       plane view (coronal=1,sagittal=2,axial=3)
%   slice       slice in the selected plane of the 3D cube
%   threshold   threshold above which dose shall be displayed
%   alpha       optional argument defining the alpha value, default is 0.6.
%               To use the default when providing a custom culormap, put in
%               an empty array by [].
%   cMap        optional argument defining the colormap, default is jet
%               if you want to use the default map with the window argument
%               you can use an empty array []
%   window      optional argument defining the displayed range. default is
%               [min(doseCube(:)) max(doseCube(:))]
%
% output
%   doseHandle: handle of the plotted dose axes
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
if nargin < 7 || isempty(cMap)
    cMap = jet(64);
end
if nargin < 6 || isempty(alpha)
    alpha = 0.6;
end
if nargin < 8 || isempty(window)
    window = [min(doseCube(:)) max(doseCube(:))];
end

if plane == 1  % Coronal plane
    dose_slice = squeeze(doseCube(slice,:,:));
elseif plane == 2 % sagittal plane
    dose_slice = squeeze(doseCube(:,slice,:));
elseif plane == 3 % Axial plane
    dose_slice = squeeze(doseCube(:,:,slice));
end
axes(axesHandle)

cMapScale = size(cMap,1) - 1;
dose_rgb = ind2rgb(uint8(cMapScale*(dose_slice - window(1))/(window(2)-window(1))),cMap);

% plot dose distribution
doseHandle = image('CData',dose_rgb,'Parent',axesHandle);
set(doseHandle,'AlphaData', alpha*(double(dose_slice)>threshold));

end


