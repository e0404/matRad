function Collimation = matRad_importFieldShapes(BeamSequence, BeamSeqNames)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to import collimator shapes from a DICOM RT plan
% 
% call
%   Collimation = matRad_importFieldShapes(BeamSequence, BeamSeqNames)
%
% input
%   BeamSequence: struct containing the BeamSequence elements from the RT    
%   BeamSeqNames: cell containing the names of the elements in BeamSequence
%
% output
%   Collimation: struct with all meta information about the collimators and
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

% get collimator device types 
DeviceSeq = BeamSequence.(BeamSeqNames{1}).BeamLimitingDeviceSequence;
DeviceSeqNames = fieldnames(DeviceSeq);

% set device specific parameters
Collimation.Devices = struct;
for i = 1:length(DeviceSeqNames)
    currLimitsSeq = DeviceSeq.(DeviceSeqNames{i});
    Collimation.Devices(i).DeviceType = currLimitsSeq.RTBeamLimitingDeviceType;
    Collimation.Devices(i).NumOfLeafs = currLimitsSeq.NumberOfLeafJawPairs;
    Collimation.Devices(i).Direction = Collimation.Devices(i).DeviceType(length(Collimation.Devices(i).DeviceType)); 
    if strncmpi(Collimation.Devices(i).DeviceType,'MLC',3)
       Collimation.Devices(i).Limits = currLimitsSeq.LeafPositionBoundaries; 
    end
end

% get all field names
numOfFields = 0;
FieldSeqNames{length(BeamSeqNames)} = [];
for i = 1:length(BeamSeqNames)
    FieldSeqNames{i} = fieldnames(BeamSequence.(BeamSeqNames{i}).ControlPointSequence);
    numOfFields = numOfFields + numel(FieldSeqNames{i});
end

Collimation.numOfFields = numOfFields;

% check for following meta data in every control point sequence
% Format: 'DICOM Name Tag' 'Name in struct'; ... 
meta =  {'NominalBeamEnergy' 'Energy';'GantryAngle' 'GantryAngle';...
        'PatientSupportAngle' 'CouchAngle';'SourceToSurfaceDistance' 'SSD'};

% save information which control point sequence belongs to which beam sequence
Collimation.FieldOfBeam = NaN * ones(numOfFields,1);

counter = 0;
Collimation.Fields = struct;
% extract field information
for i = 1:length(BeamSeqNames)
    currBeamSeq = BeamSequence.(BeamSeqNames{i});
    cumWeight = 0;
    % all meta informations must be defined in the first control point sequence
    % see DICOM Doc Sec. C.8.8.14.5
    % http://dicom.nema.org/MEDICAL/Dicom/2015c/output/chtml/part03/sect_C.8.8.14.5.html
    for j = 1:length(meta)
        FieldMeta.(meta{j,2}) = currBeamSeq.ControlPointSequence.Item_1.(meta{j,1});
    end
    for j = 1:length(FieldSeqNames{i})
        counter = counter + 1;
        currFieldSeq = currBeamSeq.ControlPointSequence.(FieldSeqNames{i}{j});
        
        % get the leaf position for every device
        Collimation.Fields(counter).LeafPos{length(DeviceSeqNames),1} = [];
        for k = 1:length(DeviceSeqNames)
            currLeafPos = currFieldSeq.BeamLimitingDevicePositionSequence.(DeviceSeqNames{k}).LeafJawPositions;
            
            Collimation.Fields(counter).LeafPos{k} = NaN*ones(Collimation.Devices(k).NumOfLeafs,2);
            Collimation.Fields(counter).LeafPos{k}(:,1) = currLeafPos(1:Collimation.Devices(k).NumOfLeafs);
            Collimation.Fields(counter).LeafPos{k}(:,2) = currLeafPos(Collimation.Devices(k).NumOfLeafs+1:end); 
        end
        
        % get field meta information
        newCumWeight = currBeamSeq.ControlPointSequence.(FieldSeqNames{i}{j}).CumulativeMetersetWeight;
        Collimation.Fields(counter).Weight = newCumWeight - cumWeight;
        cumWeight = newCumWeight;
        Collimation.Fields(counter).FinalCumWeight = currBeamSeq.FinalCumulativeMetersetWeight;
        Collimation.Fields(counter).SAD = currBeamSeq.SourceAxisDistance;
        
        % other meta information is only included in all control point
        % sequences if it changes during treatment, otherwise use FieldMeta
        for k = 1:length(meta)
            if isfield(currFieldSeq,meta{k,1})
                Collimation.Fields(counter).(meta{k,2}) = currFieldSeq.(meta{k,1});
            else
                Collimation.Fields(counter).(meta{k,2}) = FieldMeta.(meta{k,2});
            end
        end
        
        
        Collimation.FieldOfBeam(counter) = i;
    end
end

% use same dimensions for field shapes as for the kernel convolution in
% photon dose calculation
convLimits = 100; % [mm]
convResolution = .5; % [mm]

% calculate field shapes from leaf positions
counter = 0;
for i = 1:length(Collimation.Fields)
    counter = counter + 1;
    shape = ones(2*convLimits/convResolution); 
    for j = 1:length(DeviceSeqNames)
        % only ASYMX/Y and MLCX/Y are yet supported due to limited example dicoms
        if strncmpi(Collimation.Devices(j).DeviceType,'ASYM',4)
            type = 1;
        elseif strncmpi(Collimation.Devices(j).DeviceType,'MLC',3)
            type = 2;
        else
            errordlg('Device type not yet supported, not all field shapes could be completely imported!');
            continue;
        end
        for k = 1:Collimation.Devices(j).NumOfLeafs
            % determine corner points of the open area
            p1 = round((Collimation.Fields(i).LeafPos{j}(k,1)+convLimits)/convResolution+1); 
            p2 = round((Collimation.Fields(i).LeafPos{j}(k,2)+convLimits)/convResolution+1);
            if type == 2
                p3 = round((Collimation.Devices(j).Limits(k)+convLimits)/convResolution+1);
                p4 = round((Collimation.Devices(j).Limits(k+1)+convLimits)/convResolution+1);
            else % for one dimensional collimation (ASMX/Y) other direction is fully open
                p3 = 1;
                p4 = 2*convLimits/convResolution;
            end

            % set elements covered by the collimator to 0
            % differentiate between x and y direction
            if (p1 > 0) && (p1 <= 2*convLimits/convResolution+1) && ...
               (p2 > 0) && (p2 <= 2*convLimits/convResolution) && ...
               (p3 > 0) && (p3 <= 2*convLimits/convResolution) && ...
               (p4 > 0) && (p4 <= 2*convLimits/convResolution+1)
                try
                    if strcmpi(Collimation.Devices(j).Direction, 'X')
                        shape(p3:(p4-1),1:(p1-1)) = 0;
                        shape(p3:(p4-1),p2:end) = 0;
                    elseif strcmpi(Collimation.Devices(j).Direction, 'Y')
                        shape(1:(p1-1),p3:(p4-1)) = 0;
                        shape(p2:end,p3:(p4-1)) = 0;
                    else
                        errordlg('Wrong collimation direction given! Not all fields could be imported.');
                        continue;
                    end
                catch
                    warning('Error in setting field shapes. Only open field sizes up to 200x200 mm are supported')
                    continue;
                end
            end
        end
    end
    Collimation.Fields(i).Shape = shape;
end


end
