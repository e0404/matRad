function obj = matRad_importDicomRtss(obj)
% matRad function to read the data of the selected dicomRT structure set 
% file into a matRad structure
% 
% In your object, there must be properties that contain:
%   - name of the rtss file;
%   - meta information from the dicom ct files for sanity checks.
% Optional:
%   - boolean to turn on/off visualization.
%
% Output - structure containing names, numbers, colors and coordinates 
% of the polygon segmentations.
%
% call
%   obj = matRad_importDicomRtss(obj)
%
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

matRad_cfg.dispInfo('\nReading structures...');

if nargin < 3
    obj.visBool = 0;
end

matRad_checkEnvDicomRequirements(matRad_cfg.env);



% read dicom info (this includes already all data for the rtss)
if matRad_cfg.isOctave || verLessThan('matlab','9')
    structInfo = dicominfo(obj.importFiles.rtss{1});
else % apply 'UseVRHeuristic' option when available to use a to help read certain 
     % noncompliant files which switch value representation (VR) modes incorrectly
    structInfo = dicominfo(obj.importFiles.rtss{1},'UseVRHeuristic',false,'UseDictionaryVR',true);
end

% list the defined structures
try
    listOfDefStructs = fieldnames(structInfo.StructureSetROISequence);
catch
    matRad_cfg.dispError('StructureSetROISequence not defined ')
end
% list of contoured structures
try
    listOfContStructs = fieldnames(structInfo.ROIContourSequence);
catch
    matRad_cfg.dispError('ROIContourSequence not defined ')
end 

%% process structure data
numOfDefStructs  = numel(listOfDefStructs);
numOfContStructs = numel(listOfContStructs);

for i = 1:numOfContStructs % loop over every structure   

    % find the correct name
    for j = 1:numOfDefStructs
        if structInfo.ROIContourSequence.(listOfContStructs{i}).ReferencedROINumber ...
                == structInfo.StructureSetROISequence.(listOfDefStructs{j}).ROINumber
            break;
        end
    end    
    obj.importRtss.structures(i).structName   = regexprep(...  % replace nonregular characters by whitespace
        structInfo.StructureSetROISequence.(listOfDefStructs{j}).ROIName,...
        '[^a-zA-Z0-9]',' ');
                             
    obj.importRtss.structures(i).structNumber = structInfo.ROIContourSequence.(...
                                 listOfContStructs{i}).ReferencedROINumber;
    if isfield(structInfo.ROIContourSequence.(listOfContStructs{i}),'ROIDisplayColor')
        obj.importRtss.structures(i).structColor  = structInfo.ROIContourSequence.(...
                                     listOfContStructs{i}).ROIDisplayColor;  
    end

    if isfield(structInfo.ROIContourSequence.(...
                    listOfContStructs{i}), 'ContourSequence');
                if ~isempty(structInfo.ROIContourSequence.(...
                                listOfContStructs{i}).ContourSequence);
                    listOfSlices = fieldnames(structInfo.ROIContourSequence.(...
                                                listOfContStructs{i}).ContourSequence);
                else
                    matRad_cfg.dispWarning(['Contour ' obj.importRtss.structures(i).structName ' is empty'])
                    continue;
                end
    else
        matRad_cfg.dispWarning(['Contour ' obj.importRtss.structures(i).structName ' is empty'])
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
        
        % rounding to solve numerical problems with contour points not
        % being defined exactly in the same slice
        structZ = 1e-10*round(1e10*structZ);
        
        % sanity check 1
        if numel(unique(structZ)) > 1
            matRad_cfg.dispError('Detected contour points outside of single slice\n');
        end
        
        % sanity check 2
         if unique(structZ) > max(obj.ct.dicomInfo.SlicePositions) || unique(structZ) < min(obj.ct.dicomInfo.SlicePositions)
            matRad_cfg.dispWarning(['Omitting contour data for ' obj.importRtss.structures(i).structName ' at slice position ' num2str(unique(structZ)) 'mm - no ct data available.\n']);
        else
            obj.importRtss.structures(i).item(j).points = [structX, structY, structZ];
        end
            
    end
    
end

%% visualization
% show all structure points in a single plot
if obj.visBool
    figure;
    hold on
    for i = 1:numel(obj.importRtss.structures)
        plot3(obj.importRtss.structures(i).points(:,1),obj.importRtss.structures(i).points(:,2),...
            obj.importRtss.structures(i).points(:,3),'-',...
            'Color',obj.importRtss.structures(i).structColor ./ 255,'Displayname',obj.importRtss.structures(i).structName);
    end
    legend('show')
end

matRad_cfg.dispInfo('finished!\n');


end
