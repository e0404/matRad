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

    %Identifiers
    meta.SOPInstanceUID             = dicomuid;
    meta.MediaStorageSOPInstanceUID = meta.SOPInstanceUID;
    meta.SeriesInstanceUID          = dicomuid;
    meta.SeriesNumber               = 1;
    meta.InstanceNumber             = 1;

    %Remaining Meta Data
    meta.Modality = 'RTPLAN';
    meta.Manufacturer = '';
    meta.ReferringPhysicianName = obj.dicomName();
    meta.OperatorsName = obj.OperatorsName;
    meta.StationName = '';
    meta = obj.assignDefaultMetaValue(meta,'ManufacturerModelName','matRad DicomExport');

    meta.PatientName = obj.PatientName;
    meta.PatientID = obj.PatientID;
    meta.PatientBirthDate = obj.PatientBirthDate;
    meta.PatientSex = obj.PatientSex;
end


% Write Photon or ion RTPLAN
if strcmp(pln.radiationMode,'photons')
    BeamParam = 'BeamSequence';
    ControlParam = 'ControlPointSequence';
elseif strcmp(pln.radiationmode, 'protons')||strcmp(pln.radiationmode, 'carbon')||strcmp(pln.radiationmode, 'helium')
    BeamParam = 'IonBeamSequence';
    ControlParam = 'IonControlPointSequence';
else
    matRad_cfg.DispError('Not supported radiation mode of DICOM RT plan file.')
end

if strcmp(pln.radiationMode,'photons')
    %% FRaction number
    %  meta.FractionGroupSequence.Item_1.
    meta.FractionGroupSequence.Item_1.FractionGroupNumber =  1;
    meta.FractionGroupSequence.Item_1.NumberOfFractionsPlanned =  pln.numOfFractions;
    meta.FractionGroupSequence.Item_1.NumberOfBeams = pln.propStf.numOfBeams;

    %Ingnoring
    % NumberOfBrachyApplicationSetups: 0
    % ReferencedBeamSequence: [1×1 struct]
    % meta.FractionGroupSequence.Item_1.ReferencedBeamSequence.Item_1
    % BeamDoseSpecificationPoint: [3×1 double]
    %                       BeamDose: 0.3512
    %                   BeamMeterset: 116.7359
    %           ReferencedBeamNumber: 1


    % Write beam sequence info 
    if pln.propStf.numOfBeams == numel(stf)
        for i = 1: pln.propStf.numOfBeams
            %     meta.BeamSequence.(['Item_' num2str(i)]).
            meta.BeamSequence.(['Item_' num2str(i)]).RadiationType = 'PHOTON'; %upper(pln.radiationMode);
            meta.BeamSequence.(['Item_' num2str(i)]).BeamType = 'STATIC';
            meta.BeamSequence.(['Item_' num2str(i)]).PrimaryDosimeterUnit = 'MU';
            meta.BeamSequence.(['Item_' num2str(i)]).SourceAxisDistance = stf(i).SAD;
            meta.BeamSequence.(['Item_' num2str(i)]).TreatmentDeliveryType = 'TREATMENT';
            meta.BeamSequence.(['Item_' num2str(i)]).NumberOfControlPoints = 2;
            meta.BeamSequence.(['Item_' num2str(i)]).TreatmentMachineName = pln.machine;

            %       ignoring:
            %        BeamLimitingDeviceSequence: [1×1 struct]
            %                        BeamNumber: 1
            %                          BeamName: 'FB1'
            %                   BeamDescription: '1 FB1'
            %                    NumberOfWedges: 1
            %                     WedgeSequence: [1×1 struct]
            %              NumberOfCompensators: 0
            %                      NumberOfBoli: 0
            %                    NumberOfBlocks: 0
            %     FinalCumulativeMetersetWeight: 1
            %             NumberOfControlPoints: 2
            %      ReferencedPatientSetupNumber: 1



            % need a way to handle control point sequences for photons and ions
            %Beamdata
            meta.BeamSequence.(['Item_' num2str(i)]).ControlPointSequence.Item_1.ControlPointIndex = 0;
            meta.BeamSequence.(['Item_' num2str(i)]).ControlPointSequence.Item_1.NominalBeamEnergy = stf(1).ray.energy;
            meta.BeamSequence.(['Item_' num2str(i)]).ControlPointSequence.Item_1.GantryAngle = stf(i).gantryAngle;
            meta.BeamSequence.(['Item_' num2str(i)]).ControlPointSequence.Item_1.PatientSupportAngle = stf(i).couchAngle;
            meta.BeamSequence.(['Item_' num2str(i)]).ControlPointSequence.Item_1.TableTopEccentricAngle = 0;
            meta.BeamSequence.(['Item_' num2str(i)]).ControlPointSequence.Item_1.IsocenterPosition = stf(i).isoCenter';

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
