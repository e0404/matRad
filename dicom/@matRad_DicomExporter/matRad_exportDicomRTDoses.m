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
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 matRad_cfg = MatRad_Config.instance();
matRad_cfg.dispInfo('Exporting DICOM RTDose...\n');

env = matRad_getEnvironment();
isOctave = strcmp(env,'OCTAVE');

if isOctave
    matRad_cfg.dispWarning('RTDose export currently not supported by matRad running in Octave using the dicom package! Skipping...');
    return;
end

%% CT check
ct = obj.ct;
if ~any(isfield(ct,{'x','y','z'}))
    %positionOffset = transpose(ct.cubeDim ./ 2);
    positionOffset = ct.cubeDim ./ 2;
    ct.x = ct.resolution.x*[0:ct.cubeDim(2)-1] - positionOffset(2);
    ct.y = ct.resolution.y*[0:ct.cubeDim(1)-1] - positionOffset(1);
    ct.z = ct.resolution.z*[0:ct.cubeDim(3)-1] - positionOffset(3);
end

%% Meta data
storageClass = '1.2.840.10008.5.1.4.1.1.481.2';
meta.MediaStorageSOPClassUID = storageClass;
meta.SOPClassUID = storageClass;

%TransferSyntaxUID = '1.2.840.10008.1.2.1'; %Explicit VR little endian?
%meta.TransferSyntaxUID = TransferSyntaxUID;

meta.Modality = 'RTDOSE';
meta.Manufacturer = '';
meta.DoseUnits = 'GY';

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

meta.FrameOfReferenceUID = obj.FrameOfReferenceUID;
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
meta.SliceThickness = resolution.z;

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
try
    rtPlanUID = obj.rtPlanMeta.SOPInstanceUID;
    rtPlanClassID = obj.rtPlanMeta.SOPClassUID;
    meta.ReferencedRTPlanSequence.Item_1.ReferencedSOPClassUID = rtPlanClassID;
    meta.ReferencedRTPlanSequence.Item_1.ReferencedSOPInstanceUID = rtPlanUID;
catch
    rtPlanUID = '';
    rtPlanClassID = '';
end
    


if nargin < 4 || isempty(doseFieldNames)
    doseFieldNames = cell(0);
    fn = fieldnames(obj.resultGUI);
    for i = 1:numel(fn)
        if numel(size(obj.resultGUI.(fn{i}))) == 3
            doseFieldNames{end+1} = fn{i};
        end
    end
end





obj.rtDoseMetas = struct([]);
obj.rtDoseExportStatus = struct([]);

for i = 1:numel(doseFieldNames)
    doseName = doseFieldNames{i};
    
    %Now check if we export physical or RBE weighted dose, they are known
    %to dicom
    if strncmp(doseName,'physicalDose',12)
        doseType = 'PHYSICAL';
    elseif strncmp(doseName,'RBExDose',8)
        doseType = 'EFFECTIVE';
    else
        fprintf('Dose Cube ''%s'' of unknown type for DICOM. Not exported!\n',doseName);
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
        fprintf('Dose Cube ''%s'' has negative values. Not exported!\n',doseName);
        continue;
    end
    
    doseCubeFac = maxDose / double(uint16(Inf));   
    doseCube = uint16(doseCube ./ doseCubeFac);
         
    metaCube = meta;
    metaCube.DoseType = doseType;
    metaCube.DoseSummationType = deliveryType;
       
    %ID of the RTDose
    meta.SeriesInstanceUID = dicomuid;    
    metaCube.SeriesNumber = i;
    metaCube.InstanceNumber = 1;
    metaCube.DoseGridScaling = doseCubeFac;
    
    meta.SOPInstanceUID = dicomuid;
    meta.MediaStorageSOPInstanceUID = meta.SOPInstanceUID;
    
    fileName = [obj.rtDoseFilePrefix num2str(i) '_' doseName '.dcm'];
    
    env = matRad_getEnvironment();
    if strcmp(env,'OCTAVE')
        dicomwrite(doseCube,fullfile(obj.dicomDir,fileName),metaCube);
    else
        status = dicomwrite(doseCube,fullfile(obj.dicomDir,fileName),metaCube,'CreateMode','copy');%,'TransferSyntax',TransferSyntaxUID);
        if ~isempty(status)
            obj.rtDoseExportStatus = obj.addStruct2StructArray(obj.rtDoseExportStatus,status,i);
        end
    end
    
    
    obj.rtDoseMetas = obj.addStruct2StructArray(obj.rtDoseMetas,metaCube);
    
    obj.rtDoseNames{i} = doseFieldNames{i};
    
    matRad_progress(i,numel(doseFieldNames));
end
               

end
