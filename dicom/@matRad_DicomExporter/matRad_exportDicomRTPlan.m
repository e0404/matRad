function obj = matRad_exportDicomRTPlan(obj)
% matRad function to export pln to dicom.
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
    matRad_cfg.dispWarning('RTPlan export currently not supported by matRad running in Octave using the dicom package! Skipping...');
    return;
end
%%



%%

if isfield(pln,'DicomInfo')
    % if from an imported plan then use the existing dicom info
    meta = pln.DicomInfo.Meta;

else
    meta = struct([]);

    %Class UID
    ClassUID = '1.2.840.10008.5.1.4.1.1.481.5'; %RT PLAN
    meta.MediaStorageSOPClassUID = ClassUID;
    meta.SOPClassUID             = ClassUID;
    %TransferSyntaxUID = '1.2.840.10008.1.2.1'; %Explicit VR Little Endian - correct?
    %meta.TransferSyntaxUID = TransferSyntaxUID;

    %Identifiers for this object
    meta.SOPInstanceUID             = dicomuid;
    meta.MediaStorageSOPInstanceUID = meta.SOPInstanceUID;    
    meta.SeriesInstanceUID          = dicomuid;
    meta.SeriesNumber               = 1;
    meta.InstanceNumber             = 1;
    
    currDate = now;
    currDateStr = datestr(currDate,'yyyymmdd');
    currTimeStr = datestr(currDate,'HHMMSS');
    meta.InstanceCreationDate = currDateStr;
    meta.InstanceCreationTime = currTimeStr;
    
    %ID/Information about the Study
    meta.StudyInstanceUID = obj.StudyInstanceUID;
    meta.StudyID = obj.StudyID; 
    meta.StudyDate = obj.StudyDate;
    meta.StudyTime = obj.StudyTime;
      
    %Remaining Meta Data
    meta.Modality = 'RTPLAN';
    meta.Manufacturer = 'matRad';
    meta.ReferringPhysicianName = obj.dicomName();
    meta.OperatorsName = obj.OperatorsName;
    meta.StationName = '';
    meta.AccessionNumber = '';
    meta = obj.assignDefaultMetaValue(meta,'ManufacturerModelName','matRad DicomExport');

    meta.PatientName = obj.PatientName;
    meta.PatientID = obj.PatientID;
    meta.PatientBirthDate = obj.PatientBirthDate;
    meta.PatientSex = obj.PatientSex;

    meta.ApprovalStatus = 'UNAPPROVED';

    meta.RTPlanLabel = obj.rtPlanLabel;
    meta.RTPlanName  = obj.rtPlanName;
    meta.RTPlanDate = currDateStr;
    meta.RTPlanTime = currTimeStr;
    meta.RTPlanGeometry = 'PATIENT';

    meta.SoftwareVersions = matRad_version();
    meta.ManufacturerModelName = matRad_version();       
    
    %Referenced Dose
    try
        rtPlanUID = obj.rtPlanMeta.SOPInstanceUID;
        rtPlanClassID = obj.rtPlanMeta.SOPClassUID;
        meta.ReferencedRTPlanSequence.Item_1.ReferencedSOPClassUID = rtPlanClassID;
        meta.ReferencedRTPlanSequence.Item_1.ReferencedSOPInstanceUID = rtPlanUID;
    catch
        rtPlanUID = '';
        rtPlanClassID = '';
    end
end


%Write references to image, RTStruct, RTDose
%RTStruct
meta.ReferencedStructureSetSequence.Item_1.ReferencedSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3';
meta.ReferencedStructureSetSequence.Item_1.ReferencedSOPInstanceUID = obj.rtssMeta.SOPInstanceUID;

%RTDose - %TODO
%meta.DoseReferenceSquence.Item_1. ...


if pln.propStf.numOfBeams ~= numel(stf)
    matRad_cfg.error('Inconsistency in stf! number of beams not matching!');
end
   

%Sequences
%Sequences
%ToleranceTableSequence - Optional, we do not write this

 % Fraction Sequence
 %  meta.FractionGroupSequence.Item_1.
 meta.FractionGroupSequence.Item_1.FractionGroupNumber =  1;
 meta.FractionGroupSequence.Item_1.NumberOfFractionsPlanned =  pln.numOfFractions;
 meta.FractionGroupSequence.Item_1.NumberOfBeams = pln.propStf.numOfBeams;
 meta.FractionGroupSequence.Item_1.NumberOfBrachyApplicationSetups = 0;

for i = 1:pln.propStf.numOfBeams
    %PatientSetupSequence - Required
    currItemStr = sprintf('Item_%d',i);
    meta.PatientSetupSequence.(currItemStr) = struct('PatientPosition','HFS','PatientSetupNumber',i,'SetupTechniqueDescription','');
    
   
    %TODO meta.FractionGroupSequence.Item_1.ReferencedBeamSequence - %Maybe do when going through the
    %BeamSequence or afterwards?
     

    % Write Photon or ion RTPLAN
    if strcmp(pln.radiationMode,'photons')
        BeamParam = 'BeamSequence';
        ControlParam = 'ControlPointSequence';
    elseif strcmp(pln.radiationmode, 'protons') ||strcmp(pln.radiationmode, 'helium') || strcmp(pln.radiationmode, 'carbon')
        BeamParam = 'IonBeamSequence';
        ControlParam = 'IonControlPointSequence';
    else
        matRad_cfg.DispError('Not supported radiation mode of DICOM RT plan file.')
    end

    if strcmp(pln.radiationMode,'photons')
        
            seqItem.TreatmentMachineName = pln.machine;
            seqItem.PrimaryDosimeterUnit = 'MU';
            seqItem.SourceAxisDistance = stf(i).SAD;
            seqItem.BeamNumber = i;
            seqItem.BeamType = 'STATIC';
            seqItem.RadiationType = 'PHOTON'; 
            seqItem.TreatmentDeliveryType = 'TREATMENT';
            seqItem.NumberOfWedges = 0;
            seqItem.NumberOfCompensators = 0;
            seqItem.NumberOfBoli = 0;
            seqItem.NumberOfBlocks = 0;

            %TODO: seqItem.ReferencedDoseSequence (for the beam doses??)
            
            %Meterset
            % The Meterset at a given Control Point is equal to Beam 
            % Meterset (300A,0086) specified in the Referenced Beam 
            % Sequence (300C,0004) of the RT Fraction Scheme Module, 
            % multiplied by the Cumulative Meterset Weight (300A,0134) 
            % for the Control Point, divided by the Final Cumulative 
            % Meterset Weight (300A,010E). 
            % The Meterset is specified in units defined by Primary 
            % Dosimeter Unit (300A,00B3).
            
            
            
            
           
            %TODO: Get info from apertures
            seqItem.NumberOfControlPoints = 2;
            seqItem.Item_1.ControlPointIndex = 0;
            seqItem.Item_1.NominalBeamEnergy = stf(1).ray.energy;
            seqItem.Item_1.GantryAngle = stf(i).gantryAngle;
            seqItem.Item_1.PatientSupportAngle = stf(i).couchAngle;
            seqItem.Item_1.TableTopEccentricAngle = 0;
            seqItem.Item_1.IsocenterPosition = stf(i).isoCenter';
            %BBeamLimitingDeviceSequence: [1×1 struct]
            %       ignoring:
            %        
            %                        BeamNumber: 1
            %                          BeamName: 'FB1'
            %                   BeamDescription: '1 FB1'
            %                     WedgeSequence: [1×1 struct]
            %                      NumberOfBoli: 0
            %                    NumberOfBlocks: 0
            %     FinalCumulativeMetersetWeight: 1
            %             NumberOfControlPoints: 2
            %      ReferencedPatientSetupNumber: 1



            % need a way to handle control point sequences for photons and ions
            %Beamdata
            

            % Ignoring :
            %      edgePositionSequence: [1×1 struct]
            %      BeamLimitingDevicePositionSequence: [1×1 struct]
            %      GantryRotationDirection: 'NONE'
            %      BeamLimitingDeviceAngle: 270
            %      BeamLimitingDeviceRotationDirection: 'NONE'
            %      PatientSupportRotationDirection: 'NONE'
            %      TableTopEccentricRotationDirection: 'NONE'
            %      TableTopVerticalPosition: []
            %      TableTopLongitudinalPosition: []
            %      TableTopLateralPosition: []
            %      SourceToSurfaceDistance: 941.2400            WHERE is this stored ?
            %      CumulativeMetersetWeight: 0

            meta.BeamSequence.(['Item_' num2str(i)]).ControlPointSequence.Item_2.ControlPointIndex = 1;
            %   CumulativeMetersetWeight: 1
        end

    end



end
%%

filename = 'RTplan.dcm';
filepath = obj.dicomDir;
filename = fullfile(filepath,filename);

% write dicom file
env = matRad_getEnvironment();
if isOctave
    dicomwrite(int16(zeros(2)),filename,meta);
else
    obj.rtPlanExportStatus = dicomwrite([],filename,meta,'CreateMode','copy');%,'TransferSyntax',TransferSyntaxUID);
end
obj.rtPlanMetas = meta;

end
