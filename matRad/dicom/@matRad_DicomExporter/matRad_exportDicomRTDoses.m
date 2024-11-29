function obj = matRad_exportDicomRTDoses(obj)
% matRad function to exportt resultGUI to dicom RT dose. 
% Function of matRad_DicomExporter
% 
% call
%   matRad_DicomExporter.matRad_exportDicomRTDoses()
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
matRad_cfg.dispInfo('Exporting DICOM RTDose...\n');

if matRad_cfg.isOctave
    matRad_cfg.dispWarning('RTDose export currently not supported by matRad running in Octave using the dicom package! Skipping...');
    return;
end

%% CT check
ct = obj.ct;
ct = matRad_getWorldAxes(ct);


%% Meta data
storageClass = obj.rtDoseClassUID;
meta.MediaStorageSOPClassUID = storageClass;
meta.SOPClassUID = storageClass;
meta.FrameOfReferenceUID = obj.FrameOfReferenceUID;
%TransferSyntaxUID = '1.2.840.10008.1.2.1'; %Explicit VR little endian?
%meta.TransferSyntaxUID = TransferSyntaxUID;

meta.Modality = 'RTDOSE';
meta.Manufacturer = '';

%Reference
%ID of the CT
meta.StudyInstanceUID = obj.StudyInstanceUID;
meta.StudyID = obj.StudyID; 

%Dates & Times
currDate = now;
currDateStr = datestr(currDate,'yyyymmdd');
currTimeStr = datestr(currDate,'HHMMSS');
meta.InstanceCreationDate = currDateStr;
meta.InstanceCreationTime = currTimeStr;
meta.StudyDate = obj.StudyDate;
meta.StudyTime = obj.StudyTime;

meta.PositionReferenceIndicator = '';


%Remaining stuff
meta.AccessionNumber = '';
meta.StationName = '';
meta.OperatorsName = obj.OperatorsName;
meta.ReferringPhysicianName = obj.dicomName();

meta.PatientName = obj.PatientName;
meta.PatientID = obj.PatientID;
meta.PatientBirthDate = obj.PatientBirthDate;
meta.PatientSex = obj.PatientSex;

%This RTDose series
meta.SeriesInstanceUID = dicomuid;

%Now image meta
resolution = ct.resolution;
meta.PixelSpacing = [resolution.y; resolution.x];
meta.SliceThickness = num2str(resolution.z);

meta.ImagePositionPatient = [ct.x(1); ct.y(1); ct.z(1)];
meta.ImageOrientationPatient = [1;0;0;0;1;0];


%No (Re)scaling
meta.RescaleSlope = 1;
meta.RescaleIntercept = 0;
meta.RescaleType = 'US';

%RTDose is stored as multiframe image
meta.ImagesInAcquisition = ct.cubeDim(3);
meta.NumberOfFrames = ct.cubeDim(3);
meta.GridFrameOffsetVector = transpose(ct.z - ct.z(1));

%Referenced Plan
%This does currently not work well due to how Matlab creates UIDs by
%itself, we can not know the reference before it is written by the 
%RTPlanExport, which itself needs the RTDose UIDs.
%However, we need to set the ReferencedRTPlanSequence, because it is a
%conditionally required field according to the DICOM standard
try
    rtPlanUID = obj.rtPlanMeta.SOPInstanceUID;   
catch
    obj.rtPlanMeta = struct();
    obj.rtPlanMeta.SOPInstanceUID = dicomuid;
    obj.rtPlanMeta.SOPClassUID = obj.rtPlanClassUID;
    rtPlanUID = obj.rtPlanMeta.SOPInstanceUID;
end
   
meta.ReferencedRTPlanSequence.Item_1.ReferencedSOPClassUID = obj.rtPlanClassUID;
meta.ReferencedRTPlanSequence.Item_1.ReferencedSOPInstanceUID = rtPlanUID;

doseFieldNames = cell(0);
fn = fieldnames(obj.resultGUI);
for i = 1:numel(fn)
    if numel(size(obj.resultGUI.(fn{i}))) == 3
        doseFieldNames{end+1} = fn{i};
    end
end

obj.rtDoseMetas = struct([]);
obj.rtDoseExportStatus = struct([]);

for i = 1:numel(doseFieldNames)
    doseName = doseFieldNames{i};    
    doseUnits = 'GY';
    
    %Now check if we export physical or RBE weighted dose, they are known
    %to dicom
    if strncmp(doseName,'physicalDose',12) || strncmp(doseName,'LET',3)
        doseType = 'PHYSICAL';
    elseif strncmp(doseName,'RBExDose',8) || strncmp(doseName,'BED',3) || strncmp(doseName,'alpha',3) || strncmp(doseName,'beta',3) || strncmp(doseName,'effect',6) || strncmp(doseName,'RBE',3)
        doseType = 'EFFECTIVE';
        if strncmp(doseName,'RBE',3)
            doseUnits = 'RELATIVE';
        end
    else
        matRad_cfg.dispInfo('Dose Cube ''%s'' of unknown type for DICOM. Not exported!\n',doseName);
        continue;
    end
    
    %Now check if we export a single beam
    if ~isempty(regexp(doseName,'_beam(\d+)','once'))
        %We export a single beam fraction
        deliveryType = 'BEAM_SESSION';
    else
        deliveryType = 'FRACTION';
    end
    
    
    doseCube = zeros(ct.cubeDim(1),ct.cubeDim(2),1,ct.cubeDim(3));
    doseCube(:,:,1,:) = obj.resultGUI.(doseFieldNames{i});
    
    minDose = min(doseCube(:));
    maxDose = max(doseCube(:));
    
    if minDose < 0
        matRad_cfg.dispInfo('Dose Cube ''%s'' has negative values. Not exported!\n',doseName);
        continue;
    end
    
    doseCubeFac = maxDose / double(uint16(Inf));   
    doseCube = uint16(doseCube ./ doseCubeFac);
         
    metaCube = meta;
    metaCube.DoseType = doseType;
    metaCube.DoseSummationType = deliveryType;
    metaCube.DoseComment = doseName;
    metaCube.DoseUnits = doseUnits;
       
    %ID of the RTDose
    metaCube.SeriesInstanceUID = dicomuid;    
    metaCube.SeriesNumber = i;
    metaCube.InstanceNumber = 1;
    metaCube.DoseGridScaling = doseCubeFac;
    
    metaCube.SOPInstanceUID = dicomuid;
    metaCube.MediaStorageSOPInstanceUID = metaCube.SOPInstanceUID;
    
    fileName = [obj.rtDoseFilePrefix num2str(i) '_' doseName '.dcm'];    
    fileName = fullfile(obj.dicomDir,fileName);

    if matRad_cfg.isOctave
        dicomwrite(doseCube,fileName,metaCube);
    else
        status = dicomwrite(doseCube,fileName,metaCube,'CreateMode','copy');%,'TransferSyntax',TransferSyntaxUID);
        if ~isempty(status)
            obj.rtDoseExportStatus = obj.addStruct2StructArray(obj.rtDoseExportStatus,status,i);
        end
    end

    %We need to get the info of the file just written because of Matlab's
    %hardcoded way of generating InstanceUIDs during writing
    tmpInfo = dicominfo(fileName);
    metaCube.SOPInstanceUID              = tmpInfo.SOPInstanceUID;
    metaCube.MediaStorageSOPInstanceUID  = tmpInfo.MediaStorageSOPInstanceUID;   
    
    
    obj.rtDoseMetas = obj.addStruct2StructArray(obj.rtDoseMetas,metaCube);
    
    obj.rtDoseNames{i} = doseFieldNames{i};
    
    % Show progress
    if matRad_cfg.logLevel > 2
        matRad_progress(i,numel(doseFieldNames));
    end
    
end
               

end
