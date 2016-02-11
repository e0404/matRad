function structures = matRad_importDicomRtss(filename,dicomInfo,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to read the data of the selected dicomRT structure set 
% file into a matlab struct
% 
% call
%   structures = matRad_importDicomRtss(filename,dicomInfo,visBool)
%
% input
%   filename:       name of the rtss file
%   dicomInfo:      meta information from the dicom ct files for sanity
%                   checks
%   visBool:        optional: turn on/off visualization
%
% output
%   structures:     struct containing names, numbers, colors, and
%                   coordinates of the polygon segmentations
%
% References
%   -
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

fprintf('\nReading structures...');

if nargin < 3
    visBool = 0;
end

%% get info
structInfo = dicominfo(filename);

% list the defined structures
listOfDefStructs = fieldnames(structInfo.StructureSetROISequence);
% list of contoured structures
listOfContStructs = fieldnames(structInfo.ROIContourSequence);

%% process structure data
% numOfDefStructs = numel(listOfDefStructs);
numOfContStructs = numel(listOfContStructs);

for i = 1:numOfContStructs % loop over every structure   

% Is This enough?
    structures(i).structName   = structInfo.StructureSetROISequence.(...
                                 listOfDefStructs{i}).ROIName;                               
    structures(i).structNumber = structInfo.ROIContourSequence.(...
                                 listOfContStructs{i}).ReferencedROINumber;
    structures(i).structColor  = structInfo.ROIContourSequence.(...
                                 listOfContStructs{i}).ROIDisplayColor;  
                             
    if isfield(structInfo.ROIContourSequence.(...
                    listOfContStructs{i}), 'ContourSequence');
                if ~isempty(structInfo.ROIContourSequence.(...
                                listOfContStructs{i}).ContourSequence);
                    listOfSlices = fieldnames(structInfo.ROIContourSequence.(...
                                                listOfContStructs{i}).ContourSequence);
                else
                    warning(['Contour ' structures(i).structName ' is empty'])
                    continue;
                end
    else
        warning(['Contour ' structures(i).structName ' is empty'])
        continue;
    end
    
    for j = 1:numel(listOfSlices)
        structSlice = structInfo.ROIContourSequence.(...
                listOfContStructs{i}).ContourSequence.(listOfSlices{j});
        
        if strcmpi(structSlice.ContourGeometricType, 'POINT')
            continue;
        end
        
        % store the z-coordinate of this structure slice
        structX = structSlice.ContourData([1:3:end 1]);
        structY = structSlice.ContourData([2:3:end 2]);
        structZ = structSlice.ContourData([3:3:end 3]);
        
        % sanity check 1
        if numel(unique(structZ)) > 1
            error('Detected contour points outside of single slice\n');
        end
        
        % sanity check 2
        if unique(structZ) > max(dicomInfo.SlicePositions) || unique(structZ) < min(dicomInfo.SlicePositions)
            warning(['Omitting contour data for ' structures(i).structName ' at slice position ' num2str(unique(structZ)) 'mm - no ct data available.\n']);
        else
            structures(i).item(j).points = [structX, structY, structZ];
        end
            
    end
    
end

%% visualization
% show all structure points in a single plot
if visBool
    figure;
    hold on
    for i = 1:numel(structures)
        plot3(structures(i).points(:,1),structures(i).points(:,2),...
            structures(i).points(:,3),'-',...
            'Color',structures(i).structColor ./ 255,'Displayname',structures(i).structName);
    end
    legend('show')
end

fprintf('finished!\n');


end