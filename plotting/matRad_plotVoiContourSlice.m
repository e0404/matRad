function [voiContourHandles] = matRad_plotVoiContourSlice(axesHandle,cst,ct,ctIndex,selection,plane,slice,cMap,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function that plots the contours of the segmentations given in cst
%
% call
%   voiContourHandles = matRad_plotVoiContourSlice(axesHandle,cst,plane,slice,cMap)
%
% input
%   axesHandle          handle to axes the slice should be displayed in
%   cst                 matRad cst cell array
%   ct                  matRad ct structure
%   ctIndex             index of the ct cube
%   selection           logicals defining the current selection of contours
%                       that should be plotted. Can be set to [] to plot
%                       all non-ignored contours.
%   plane               plane view (coronal=1,sagittal=2,axial=3)
%   slice               slice in the selected plane of the 3D cube
%   cMap                optional argument defining the colormap, default is
%                       colorcube
%   varargin            Additional Matlab Line-Property/value pairs
%
% output
%   voiContourHandles:  handles of the plotted contours
%   visibleOnSlice:     logicals defining if the contour is actually
%                       visible on the current slice
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

if ~isdeployed
    addpath('tools')
end

[env, ~] = matRad_getEnvironment();

% overwrite colormap
if exist('cMap', 'var') && ~isempty(cMap)
    cMapScale = size(cMap,1)-1;
    %determine colors
    colors = cMap(round(linspace(1,cMapScale,size(cst,1))),:);
  else
    for i = 1:size(cst,1)
      colors(i,:) = cst{i,5}.visibleColor;
    end
end

if isempty(selection) || numel(selection) ~= size(cst,1)
    selection = true(size(cst,1),1);
end

voiContourHandles = cell(0);
switch env
    case 'MATLAB'
        voiContourHandles = gobjects(0);
    case 'OCTAVE'
        voiContourHandles = [];
end


for s = 1:size(cst,1)
    
    if ~strcmp(cst{s,3},'IGNORED') && selection(s)
        %Check for precalculated contours
        C =[];
        if size(cst,2) >= 7 && ~isempty(cst{s,7})
            C = cst{s,7}{slice,plane};
        else
            %If we do not have precomputed contours available, then compute them
            mask = zeros(size(ct{ctIndex}));
            mask(cst{s,4}{ctIndex}) = 1;
            
            if plane == 1 && any(any(mask(slice,:,:) > 0))
                C = contourc(squeeze(mask(slice,:,:)),.5*[1 1]);
            elseif plane == 2 && any(any(mask(:,slice,:) > 0))
                C = contourc(squeeze(mask(:,slice,:)),.5*[1 1]);
            elseif plane == 3 && any(any(mask(:,:,slice) > 0))
                C = contourc(squeeze(mask(:,:,slice)),.5*[1 1]);
            end  
        end
        
        % plot precalculated contourc data
        
         switch env
            case 'MATLAB'
                tmpLineHandle = gobjects(0);
            case 'OCTAVE'
                tmpLineHandle = [];
         end

         if any(C(:))
            lower = 1; % lower marks the beginning of a section
            while lower-1 ~= size(C,2)
                hold on
                steps = C(2,lower); % number of elements of current line section
                tmpLineHandle(end+1) = line(C(1,lower+1:lower+steps),...
                    C(2,lower+1:lower+steps),...
                    'Color',colors(s,:),'Parent',axesHandle,varargin{:});

                lower = lower+steps+1;
            end
            voiContourHandles{end+1} = tmpLineHandle;
         else
            % create empty line object
            voiContourHandles{end+1} = {};
         end     

               
    end
end


end
