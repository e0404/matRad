function doseHandle = matRad_plotDoseSlice(axesHandle,doseCube,plane,slice,threshold,alpha,cMap)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function that generates the plot for the CT in the GUI. The
% function can also be used in personal matlab figures by passing the
% corresponding axes handle
%
% call
%   doseHandle = matRad_plotDoseSlice(axesHandle,ctCube,plane,slice,cMap)
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
%
% output
%   doseHandle:   handle of the plotted dose axes
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

if plane == 1  % Coronal plane
    dose_slice = squeeze(doseCube(slice,:,:));
elseif plane == 2 % sagittal plane
    dose_slice = squeeze(doseCube(:,slice,:));
elseif plane == 3 % Axial plane
    dose_slice = squeeze(doseCube(:,:,slice));
end
axes(axesHandle)

cMapScale = size(cMap,1) - 1;
dose_rgb = ind2rgb(uint8(cMapScale*dose_slice),cMap);

% plot dose distribution
doseHandle = imagesc('CData',dose_rgb,'Parent',axesHandle);
set(doseHandle,'AlphaData', alpha*(double(dose_slice)>threshold));

end


