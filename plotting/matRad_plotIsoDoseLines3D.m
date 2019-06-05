function isoLineHandles = matRad_plotIsoDoseLines3D(axesHandle,ct,doseCube,isoContours,isoLevels,plane,slice,cMap,window,varargin)
% matRad function that plots isolines in 3D, by precomputed contourc data 
% computed by matRad_computeIsoDoseContours or manually by calling 
% contourslice itself
%
% call
%   isoLineHandles = matRad_plotIsoDoseLines3D(axesHandle,ct,doseCube,isoContours,isoLevels,plotLabels,plane,slice,cMap,window)
%
% input
%   axesHandle  handle to axes the slice should be displayed in
%   ct          matRad ct struct which contains resolution
%   doseCube    3D array of the corresponding dose cube
%   isoContours precomputed isodose contours in a cell array {maxDim,3}
%               if the parameter is empty, contours will be plotted the
%               slow way with MATLABs contour function
%   isoLevels   the levels of the isodose (same units as doseCube)
%   plane       plane view (coronal=1,sagittal=2,axial=3)
%   slice       slice in the selected plane of the 3D cube
%   cMap        optional argument defining the colormap, default is jet
%               if you want to use the default map with the window argument
%               you can use an empty array []
%   window      optional argument defining the displayed range. default is
%               [min(doseCube(:)) max(doseCube(:))]
%   varargin    Additional Matlab Line-Property/value pairs
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
if nargin < 7 || isempty(cMap)
    cMap = jet(64);
end
if nargin < 8 || isempty(window)
    window = [min(doseCube(:)) max(doseCube(:))];
end

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

slices = {[],[],[]};

coords{1} = ct.resolution.x * double(1:ct.cubeDim(2));
coords{2} = ct.resolution.y * double(1:ct.cubeDim(1));
coords{3} = ct.resolution.z * double(1:ct.cubeDim(3));

%slice spacing
spacing = 5;
sliceIndices = [1:spacing:ct.cubeDim(plane)];
if isempty(sliceIndices(sliceIndices==slice))
    sliceIndices(end+1) = slice;
    sort(sliceIndices);
end
slices{plane} = coords{plane}(sliceIndices);

%% Plotting
%Check if precomputed contours where passed, if not, calculate it on the
%fly
if isempty(isoContours)
    
    [xMesh,yMesh,zMesh] = meshgrid(coords{1},coords{2},coords{3});
    isoLineHandles = contourslice(axesHandle,xMesh,yMesh,zMesh,doseCube,slices{[1 2 3]},isoLevels);
else  
    axes(axesHandle);
    hold on;
    
    for s = 1:numel(sliceIndices)
        currSlice = sliceIndices(s);
        currSlicePlaneCoords = slices{plane}(s);
        
        %opacity of the isolines. Will be fully opaque if we are on the
        %requested (and shown) slice
        if currSlice == slice
            opacity = 1;
        else
            opacity = 0.4;
        end
        
        %Check if there is a contour in the plane
        if any(isoContours{currSlice,plane}(:))
            % plot precalculated contourc data            
            lower = 1; % lower marks the beginning of a section
            while lower-1 ~= size(isoContours{currSlice,plane},2)
                steps = isoContours{currSlice,plane}(2,lower); % number of elements of current line section
                if numel(unique(isoLevels)) > 1
                    color = colors(isoLevels(:) == isoContours{currSlice,plane}(1,lower),:);
                else
                    color = unique(colors,'rows');
                end              
                
                % Align 2D Contours in3D
                isoLine2Dx = isoContours{currSlice,plane}(1,lower+1:lower+steps);
                isoLine2Dy = isoContours{currSlice,plane}(2,lower+1:lower+steps);
                if plane == 2
                    isoLine3Dx = currSlicePlaneCoords*ones(1,numel(isoLine2Dx));
                    isoLine3Dz = interp1(double(1:ct.cubeDim(3)),coords{3},isoLine2Dx);
                    isoLine3Dy = interp1(double(1:ct.cubeDim(1)),coords{1},isoLine2Dy);
                elseif plane == 1
                    isoLine3Dy = currSlicePlaneCoords*ones(1,numel(isoLine2Dx));
                    isoLine3Dx = interp1(double(1:ct.cubeDim(2)),coords{2},isoLine2Dy);
                    isoLine3Dz = interp1(double(1:ct.cubeDim(3)),coords{3},isoLine2Dx);
                elseif plane == 3
                    isoLine3Dz = currSlicePlaneCoords*ones(1,numel(isoLine2Dx));
                    isoLine3Dx = interp1(double(1:ct.cubeDim(2)),coords{2},isoLine2Dx);
                    isoLine3Dy = interp1(double(1:ct.cubeDim(1)),coords{1},isoLine2Dy);
                else 
                    continue;
                end
               
                %We render the isodose lines transparent by adding a fourth color value (undocumented)
                if verLessThan('matlab','8.5')
                    isoLineHandles(end+1) = line(isoLine3Dx,isoLine3Dy,isoLine3Dz,'Color',[color        ],'Parent',axesHandle,varargin{:});
                else
                    isoLineHandles(end+1) = line(isoLine3Dx,isoLine3Dy,isoLine3Dz,'Color',[color opacity],'Parent',axesHandle,varargin{:});
                end
                %We do not plot labels
                %{
                if plotLabels
                    text(isoLine3Dx,isoLine3Dy,isoLine3Dz,num2str(isoContours{currSlice,plane}(1,lower)),'Parent',axesHandle);
                end
                %}
                lower = lower+steps+1;
                
            end
        end
    end
    
    hold off;
end

end
