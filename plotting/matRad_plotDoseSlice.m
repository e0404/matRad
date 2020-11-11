function [doseHandle,cMap,window] = matRad_plotDoseSlice(axesHandle,doseCube,plane,slice,threshold,alpha,cMap,window)
% matRad function that generates a dose plot of a selected slice 
% The function can also be used in personal matlab figures by passing the
% corresponding axes handle.
%
% call
%   [doseHandle,cMap,window] = matRad_plotDoseSlice(axesHandle, doseCube,plane,slice,threshold)
%   [doseHandle,cMap,window] = matRad_plotDoseSlice(axesHandle, doseCube,plane,slice,threshold,alpha)
%   [doseHandle,cMap,window] = matRad_plotDoseSlice(axesHandle, doseCube,plane,slice,threshold,cMap)
%   [doseHandle,cMap,window] = matRad_plotDoseSlice(axesHandle, doseCube,plane,slice,threshold,window)
%   [doseHandle,cMap,window] = matRad_plotDoseSlice(axesHandle, doseCube,plane,slice,threshold,alpha,cMap)
%   [doseHandle,cMap,window] = matRad_plotDoseSlice(axesHandle, doseCube,plane,slice,threshold,alpha,window)
%   [doseHandle,cMap,window] = matRad_plotDoseSlice(axesHandle, doseCube,plane,slice,threshold,cMap,window)
%   [doseHandle,cMap,window] = matRad_plotDoseSlice(axesHandle, doseCube,plane,slice,threshold,alpha,cMap,window)
%
% input
%   axesHandle  handle to axes the slice should be displayed in
%   doseCube    3D array of the dose to select the slice from
%   plane       plane view (coronal=1,sagittal=2,axial=3)
%   slice       slice in the selected plane of the 3D cube
%   threshold   threshold above which the dose shall be displayed
%               for negative values (i.e. difference maps), also the values
%               smaller than the negative threshold will be displayed
%               if empty, no threshold will be applied
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
if nargin < 7 || isempty(cMap)
    cMap = jet(64);
end
if nargin < 6 || isempty(alpha)
    alpha = 0.6;
end
if nargin < 8 || isempty(window)
    window = [min(doseCube(:)) max(doseCube(:))];
end

maxDose = max(doseCube(:));

cMapScale = size(cMap,1) - 1;

if plane == 1  % Coronal plane
    dose_slice = squeeze(doseCube(slice,:,:));
elseif plane == 2 % sagittal plane
    dose_slice = squeeze(doseCube(:,slice,:));
elseif plane == 3 % Axial plane
    dose_slice = squeeze(doseCube(:,:,slice));
else
	matRad_cfg.dispError('Invalid plane ''%d'' selected for visualization!',plane);
end

if ~isempty(threshold)
    mask = alpha * (dose_slice < window(2) & dose_slice > window(1) & dose_slice > threshold*maxDose);
else
    mask = alpha * (dose_slice < window(2) & dose_slice > window(1));
end

dose_slice = uint8(cMapScale*(dose_slice - window(1))/(window(2)-window(1)));

%This circumenvents a bug in Octave when the index in the image hase the maximum value of uint8
if matRad_cfg.isOctave
	dose_slice(dose_slice == 255) = 254;
end

dose_rgb = ind2rgb(dose_slice,cMap);

% plot dose distribution
doseHandle = image('CData',dose_rgb,'Parent',axesHandle);

% alphadata for image objects is not yet implemented in Octave
set(doseHandle,'AlphaData',mask);

end


