function collimation = matRad_importFieldShapes(beamSequence, fractionSequence)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to import collimator shapes from a DICOM RT plan
% 
% call
%   collimation = matRad_importFieldShapes(beamSequence, fractionSequence)
%
% input
%   beamSequence: struct containing the beamSequence elements from the RT plan    
%   fractionSequence: struct containing the fractionGroupSequence elements from the RT plan    
%
% output
%   collimation: struct with all meta information about the collimators and
%   all field shape matrices 
%
% References
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
counter = 0;
maximumExtent = 0;
tmpCollimation.Fields = struct;

% check for following meta data in every control point sequence
% Format: 'DICOM Name Tag' 'Name in struct'; ... 
meta =  {'NominalBeamEnergy' 'Energy';'GantryAngle' 'GantryAngle';...
        'PatientSupportAngle' 'CouchAngle';'SourceToSurfaceDistance' 'SSD'};

% extract field information
beamSeqNames = fields(beamSequence);
for i = 1:length(beamSeqNames)
    
    currBeamSeq = beamSequence.(beamSeqNames{i});
    cumWeight = 0;
    
    % get total MU applied by beam i
    tmpCollimation.beamMeterset(i) = fractionSequence.ReferencedBeamSequence.(beamSeqNames{i}).BeamMeterset;
    
    % get collimator device types 
    currDeviceSeq = beamSequence.(beamSeqNames{i}).BeamLimitingDeviceSequence;
    currDeviceSeqNames = fieldnames(currDeviceSeq);
    
    % set device specific parameters
    device = struct;
    for j = 1:length(currDeviceSeqNames)
        currLimitsSeq = currDeviceSeq.(currDeviceSeqNames{j});
        device(j).DeviceType = currLimitsSeq.RTBeamLimitingDeviceType;
        device(j).NumOfLeafs = currLimitsSeq.NumberOfLeafJawPairs;
        device(j).Direction = device(j).DeviceType(end); 
        if strncmpi(device(j).DeviceType,'MLC',3)
           device(j).Limits = currLimitsSeq.LeafPositionBoundaries; 
        end
    end
    tmpCollimation.Devices{i} = device;
    
    currControlPointSeqNames = fieldnames(currBeamSeq.ControlPointSequence);
    % all meta informations must be defined in the first control point sequence
    % see DICOM Doc Sec. C.8.8.14.5
    % http://dicom.nema.org/MEDICAL/Dicom/2015c/output/chtml/part03/sect_C.8.8.14.5.html
    for j = 1:length(meta)
        try
            FieldMeta.(meta{j,2}) = currBeamSeq.ControlPointSequence.(currControlPointSeqNames{1}).(meta{j,1});
        catch
            warning(['Field ' meta{j,1} ' not found on beam sequence ' beamSeqNames{i} ...
                     '. No field shape import performed!']);
            return;
        end
    end
    
    for j = 1:length(currControlPointSeqNames)
       counter = counter + 1;
       currControlPointElement = currBeamSeq.ControlPointSequence.(currControlPointSeqNames{j});      
        
       if isfield(currControlPointElement, 'BeamLimitingDevicePositionSequence')
           % get the leaf position for every device
           tmpCollimation.Fields(counter).LeafPos{length(currDeviceSeqNames),1} = [];

           % beam limiting device position sequence has to be defined on
           % the first control point and has to be defined on following
           % points only if it changes -> default initilation if counter > 1
           if counter > 1
               for k = 1:length(currDeviceSeqNames)
                   tmpCollimation.Fields(counter).LeafPos{k} = tmpCollimation.Fields(counter-1).LeafPos{k};
               end
           end

           for k = 1:length(currDeviceSeqNames)

               if isfield(currControlPointElement.BeamLimitingDevicePositionSequence,currDeviceSeqNames{k})
                   currLeafPos = currControlPointElement.BeamLimitingDevicePositionSequence.(currDeviceSeqNames{k}).LeafJawPositions;          

                   deviceIx = find(strcmp({device(:).DeviceType}, ...
                       currControlPointElement.BeamLimitingDevicePositionSequence.(currDeviceSeqNames{k}).RTBeamLimitingDeviceType));

                   if (length(currLeafPos) ~= 2 * device(deviceIx).NumOfLeafs)
                       warning(['Number of leafs/jaws does not match given number of leaf/jaw positions in control point sequence ' ...
                                currControlPointSeqNames{j} ' on beam sequence ' beamSeqNames{i} ' for device ' ...
                                device(deviceIx).DeviceType '. No field shape import performed!']);
                       return;
                   end
                   
                   % set left and right leaf positions
                   tmpCollimation.Fields(counter).LeafPos{deviceIx}(:,1) = currLeafPos(1:device(deviceIx).NumOfLeafs);
                   tmpCollimation.Fields(counter).LeafPos{deviceIx}(:,2) = currLeafPos(device(deviceIx).NumOfLeafs+1:end);
                   % find the total maximum extent of one beam (in any direction) 
                   maximumExtent = max(maximumExtent, max(abs(currLeafPos))); % check opening direction
                   % check direction perpendicular to the openening for MLC
                   if strncmpi(device(k).DeviceType,'MLC',3)
                       maximumExtent = max(maximumExtent,max(abs(device(deviceIx).Limits)));
                   end
               end
           end
       else
           tmpCollimation.Fields(counter) = tmpCollimation.Fields(counter - 1);
       end
       
       % get field meta information
       if isfield(currControlPointElement, 'CumulativeMetersetWeight')      
           newCumWeight = currControlPointElement.CumulativeMetersetWeight;
           tmpCollimation.Fields(counter).Weight = (newCumWeight - cumWeight) / ...
                                    currBeamSeq.FinalCumulativeMetersetWeight * ...
                                    tmpCollimation.beamMeterset(i)/100;
           cumWeight = newCumWeight;
       else
           warning(['No CumulativeMetersetWeight found in control point sequence ' currControlPointSeqNames{j} ...
                    ' on beam ' beamSeqNames{i} '. No field shape import performed!']);
           return;
       end
       tmpCollimation.Fields(counter).SAD = currBeamSeq.SourceAxisDistance;
        
       % other meta information is only included in all control point
       % sequences if it changes during treatment, otherwise use FieldMeta
       for k = 1:length(meta)
           if isfield(currControlPointElement,meta{k,1})
               tmpCollimation.Fields(counter).(meta{k,2}) = currControlPointElement.(meta{k,1});
           else
               tmpCollimation.Fields(counter).(meta{k,2}) = FieldMeta.(meta{k,2});
           end
       end
       % save information which control point sequence belongs to which beam sequence       
       tmpCollimation.Fields(counter).BeamIndex = i;
    end
end
tmpCollimation.numOfFields = counter;

% field import works only if the leaf width is a multiple of the conv
% resolution
convResolution = .5; % [mm]
tmpCollimation.convResolution = convResolution;

% get temporary shape limits to calculate the shapes
shapeLimit = ceil(maximumExtent / convResolution);

% calculate field shapes from leaf positions
maximumVoxelExtent = 0;
[X,Y] = meshgrid(-shapeLimit:shapeLimit-1);
for i = 1:length(tmpCollimation.Fields)
    shape = ones(2*shapeLimit); 
    beamIndex = tmpCollimation.Fields(i).BeamIndex;
    for j = 1:length(tmpCollimation.Devices{beamIndex})
        % check for ASYM and SYM jaws == type 1
        if strncmpi(tmpCollimation.Devices{beamIndex}(j).DeviceType,'ASYM',4)
            type = 1;
        elseif (strcmpi(tmpCollimation.Devices{beamIndex}(j).DeviceType,'X') || ...
                strcmpi(tmpCollimation.Devices{beamIndex}(j).DeviceType,'Y'))
            type = 1;
        % MLC == type 2
        elseif strncmpi(tmpCollimation.Devices{beamIndex}(j).DeviceType,'MLC',3)
            type = 2;
        else
            warning(['Device type ' tmpCollimation.Devices{beamIndex}(j).DeviceType ...
                    ' not supported. Field shapes could not be imported!']);
            return;
        end
        for k = 1:tmpCollimation.Devices{beamIndex}(j).NumOfLeafs
            % determine corner points of the open area
            p1 = ceil(tmpCollimation.Fields(i).LeafPos{j}(k,1)/convResolution)+shapeLimit;
            p2 = ceil(tmpCollimation.Fields(i).LeafPos{j}(k,2)/convResolution)+shapeLimit+1;
            if type == 2
                p3 = ceil(tmpCollimation.Devices{beamIndex}(j).Limits(k)/convResolution)+shapeLimit+1;
                p4 = ceil(tmpCollimation.Devices{beamIndex}(j).Limits(k+1)/convResolution)+shapeLimit;
            else % for one dimensional collimation (ASMX/Y) other direction is fully open
                p3 = 1;
                p4 = 2*shapeLimit;
            end

            % set elements covered by the collimator to 0
            % differentiate between x and y direction
            if (p1 > 0) && (p1 <= 2*shapeLimit) && ...
               (p2 > 0) && (p2 <= 2*shapeLimit) && ...
               (p3 > 0) && (p3 <= 2*shapeLimit) && ...
               (p4 > 0) && (p4 <= 2*shapeLimit)
                try
                    if strcmpi(tmpCollimation.Devices{beamIndex}(j).Direction, 'X')
                        shape(p3:p4,1:p1) = 0;
                        shape(p3:p4,p2:end) = 0;
                    elseif strcmpi(tmpCollimation.Devices{beamIndex}(j).Direction, 'Y')
                        shape(1:p1,p3:p4) = 0;
                        shape(p2:end,p3:p4) = 0;
                    else
                        warning(['Wrong collimation direction ' tmpCollimation.Devices{beamIndex}(j).Direction ...
                                 ' given for device ' tmpCollimation.Devices{beamIndex}(j).DeviceType ...
                                 '. Fields could not be imported.']);
                        return;
                    end
                catch
                    warning('Error in setting field shapes. Field shapes could not be imported.')
                    return;
                end
            end
        end
    end
    openVoxelDistance = [X(shape == 1); Y(shape == 1)];
    maximumVoxelExtent = max(maximumVoxelExtent, max(abs(openVoxelDistance)));
    tmpCollimation.Fields(i).Shape = shape;
end

tmpCollimation.fieldWidth = 2 * maximumVoxelExtent * convResolution;
% truncate field shapes to a symmetrical field with limits maximumVoxelExtent
voxelRange = (-maximumVoxelExtent+1:maximumVoxelExtent) + shapeLimit;
for i = 1:length(tmpCollimation.Fields) 
     tmpCollimation.Fields(i).Shape = tmpCollimation.Fields(i).Shape(voxelRange, voxelRange);
end
collimation = tmpCollimation;
end
