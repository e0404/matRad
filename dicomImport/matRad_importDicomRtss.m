function structures = matRad_importDicomRtss(filename,dicomInfo,visBool)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to read the locations of the structures described in the
% RTSTRUCT file and return a struct containing:
% - contour-points off all defined structures (3-dim physical coordinates)
% - structure name and number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    structures(i).structName = structInfo.StructureSetROISequence.(...
                                   listOfDefStructs{i}).ROIName;
    structures(i).structNumber = structInfo.ROIContourSequence.(...
                                 listOfContStructs{i}).ReferencedROINumber;
    structures(i).structColor = structInfo.ROIContourSequence.(...
                                 listOfContStructs{i}).ROIDisplayColor;                         

    listOfSlices = fieldnames(structInfo.ROIContourSequence.(...
                                   listOfContStructs{i}).ContourSequence);
    
    % getting data of all structure slices
    structures(i).points = [];
    
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
            structures(i).points = vertcat(structures(i).points, [structX, structY, structZ]);
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