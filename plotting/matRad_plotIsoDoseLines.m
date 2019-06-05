function isoLineHandles = matRad_plotIsoDoseLines(axesHandle,doseCube,isoContours,isoLevels,plotLabels,plane,slice,cMap,window,varargin)
% matRad function that plots isolines, by precomputed contourc data 
% computed by matRad_computeIsoDoseContours or manually by calling contourc
% itself
%
% call
%   isoLineHandles = matRad_plotIsoDoseLines(axesHandle,doseCube,isoContours,isoLevels,plotLabels,plane,slice,cMap,window)
%
% input
%   axesHandle  handle to axes the slice should be displayed in
%   doseCube    3D array of the corresponding dose cube
%   isoContours precomputed isodose contours in a cell array {maxDim,3}
%               if the parameter is empty, contours will be plotted the
%               slow way with MATLABs contour function
%   isoLevels   the levels of the isodose (same units as doseCube)
%   plotLabels  if set to true labels will be added to the contours
%   plane       plane view (coronal=1,sagittal=2,axial=3)
%   slice       slice in the selected plane of the 3D cube
%   cMap        optional argument defining the colormap, default is jet
%               if you want to use the default map with the window argument
%               you can use an empty array []
%   window      optional argument defining the displayed range. default is
%               [min(doseCube(:)) max(doseCube(:))]
%   varargin    Additional MATLAB Line-Property/Value-Pairs etc.
%
% output
%   isoLineHandles: handle to the plotted isolines
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

[env, ~] = matRad_getEnvironment();

%% manage optional arguments
%Use default colormap?
if nargin < 8 || isempty(cMap)
    cMap = jet(64);
end
if nargin < 9 || isempty(window)
    window = [min(doseCube(:)) max(doseCube(:))];
end

%Check if precomputed contours where passed, if not, calculate it on the
%fly
if isempty(isoContours)
    if plane == 1
        C = contourc(doseCube(slice,:,:),isoLevels);
    elseif plane == 2
        C = contourc(doseCube(:,slice,:),isoLevels);
    elseif plane == 3
        C = contourc(doseCube(:,:,slice),isoLevels);
    end    
    isoContours{slice,plane} = C;
end

%% Plotting
cMapScale = size(cMap,1) - 1;
isoColorLevel = (isoLevels - window(1))./(window(2)-window(1));
isoColorLevel(isoColorLevel < 0) = 0;
isoColorLevel(isoColorLevel > 1) = 0;
colors = squeeze(ind2rgb(uint8(cMapScale*isoColorLevel),cMap));

switch env
    case 'MATLAB'
        isoLineHandles = gobjects(0);
    case 'OCTAVE'
        isoLineHandles = [];
end

axes(axesHandle);
hold on;

%Check if there is a contour in the plane
if any(isoContours{slice,plane}(:))
    % plot precalculated contourc data
    
    lower = 1; % lower marks the beginning of a section
    while lower-1 ~= size(isoContours{slice,plane},2)
        steps = isoContours{slice,plane}(2,lower); % number of elements of current line section
        if numel(unique(isoLevels)) > 1
            color = colors(isoLevels(:) == isoContours{slice,plane}(1,lower),:);
        else
            color = unique(colors,'rows'); 
        end
        isoLineHandles(end+1) = line(isoContours{slice,plane}(1,lower+1:lower+steps),...
            isoContours{slice,plane}(2,lower+1:lower+steps),...
            'Color',color,'Parent',axesHandle,varargin{:});
        if plotLabels
            text(isoContours{slice,plane}(1,lower+1),...
                isoContours{slice,plane}(2,lower+1),...
                num2str(isoContours{slice,plane}(1,lower)),'Parent',axesHandle)
        end
        lower = lower+steps+1;
        
    end
end

hold off;

end
