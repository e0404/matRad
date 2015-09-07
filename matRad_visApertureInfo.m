function matRad_visApertureInfo(apertureInfo,mode)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to visualize aperture shapes stored as struct
%
% call
%   matRad_visApertureInfo(apertureInfo,mode)
%
% input
%   apertureInfo: aperture weight and shape info struct
%   mode:         switch to display leaf numbers ('leafNum') or physical
%                 coordinates of the leaves ('physical')
%
% output
%   -
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 % set default mode to physical
    mode = 'leafNum'; % options: 'physical','leafNum'
end

% global parameters
numOfBeams = size(apertureInfo.beam,2);
bixelWidth = apertureInfo.bixelWidth;

% custom colormap
color = [0.2:0.01:0.8; 0.2:0.01:0.8; 0.2:0.01:0.8]';
color = flipud(color);
color(:,3) = 0;
color(:,2) = 0;

% loop over all beams
for i=1:numOfBeams

    % open new figure for every beam
    figure

    % get the MLC dimensions for this beam
    minX = apertureInfo.beam(i).MLCWindow(1);
    maxX = apertureInfo.beam(i).MLCWindow(2);    
    
    %get maximum weight
    wMax = max([apertureInfo.beam(i).shape(:).weight]);
    if strcmp(mode,'leafNum')

        % get the active leaf Pairs
        % the leaf indices have to be flipped in order to fit to the order of
        % the leaf positions (1st row of leafPos is lowest row in physical
        % coordinates
        activeLeafInd = flipud(find(apertureInfo.beam(i).isActiveLeafPair));
    end
    
    subplotColumns = ceil(apertureInfo.beam(i).numOfShapes/2);
    subplotLines   = ceil(apertureInfo.beam(i).numOfShapes/subplotColumns);
    
    % loop over all shapes of the beam 
    for j = 1:apertureInfo.beam(i).numOfShapes
        
        % creating subplots
        subplot(subplotLines,subplotColumns,j)

        title(['Beam: ' num2str(i) ' Shape: ' num2str(j) ' w=' ...
                num2str(apertureInfo.beam(i).shape(j).weight,2)],...
                    'Fontsize',8)
        colorInd = ceil((apertureInfo.beam(i).shape(j).weight/wMax)*61+eps);
        set(gca,'Color',color(colorInd,:));
        
        hold on
        
        if strcmp(mode,'physical')
            % loop over all active leaf pairs
            for k = 1:apertureInfo.beam(i).numOfActiveLeafPairs
                fill([minX apertureInfo.beam(i).shape(j).leftLeafPos(k) ...
                        apertureInfo.beam(i).shape(j).leftLeafPos(k) minX],...
                        [apertureInfo.beam(i).leafPairPos(k)- bixelWidth/2 ...
                        apertureInfo.beam(i).leafPairPos(k)- bixelWidth/2 ...
                        apertureInfo.beam(i).leafPairPos(k)+ bixelWidth/2 ...
                        apertureInfo.beam(i).leafPairPos(k)+ bixelWidth/2],'b')
                 fill([apertureInfo.beam(i).shape(j).rightLeafPos(k) ...
                        maxX maxX ...
                        apertureInfo.beam(i).shape(j).rightLeafPos(k)],...
                        [apertureInfo.beam(i).leafPairPos(k)- bixelWidth/2 ...
                        apertureInfo.beam(i).leafPairPos(k)- bixelWidth/2 ...
                        apertureInfo.beam(i).leafPairPos(k)+ bixelWidth/2 ...
                        apertureInfo.beam(i).leafPairPos(k)+ bixelWidth/2],'b')                
            end
        elseif strcmp(mode,'leafNum')
            % loop over all active leaf pairs
            for k = 1:apertureInfo.beam(i).numOfActiveLeafPairs
                fill([minX apertureInfo.beam(i).shape(j).leftLeafPos(k) ...
                    apertureInfo.beam(i).shape(j).leftLeafPos(k) minX],...
                    [activeLeafInd(k) - 1/2 ...
                    activeLeafInd(k) - 1/2 ...
                    activeLeafInd(k) + 1/2 ...
                    activeLeafInd(k) + 1/2],'b') 
                fill([apertureInfo.beam(i).shape(j).rightLeafPos(k) ...
                    maxX maxX ...
                    apertureInfo.beam(i).shape(j).rightLeafPos(k)],...
                    [activeLeafInd(k) - 1/2 ...
                    activeLeafInd(k) - 1/2 ...
                    activeLeafInd(k) + 1/2 ...
                    activeLeafInd(k) + 1/2],'b')                
            end
        end
    
        axis tight
        xlabel('horiz. pos. [mm]')

        if strcmp(mode,'physical')
            ylabel('vert. pos. [mm]')
        elseif strcmp(mode,'leafNum')
            ylabel('leaf pair #')
        end
    
    end

end

end