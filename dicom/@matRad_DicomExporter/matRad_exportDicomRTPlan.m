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

%Safeguard Error
matRad_cfg.dispError('You are trying to run the RTPlan Dicom Export! This error message is a safeguard such that this export is not used unconciously and that you aware that exported RTPlans MUST NOT be used cinically. If you want to use it (at your own responsibility), uncomment this error in %s!',mfilename('fulpath'));

if isOctave
    matRad_cfg.dispWarning('RTPlan export currently not supported by matRad running in Octave using the dicom package! Skipping...');
    return;
end

%%
matRad_cfg.dispInfo('Exporting RTPlan...\n');
if isfield(obj.pln,'DicomInfo')
    % if from an imported plan then use the existing dicom info
    meta = pln.DicomInfo.Meta;

else
    meta = struct();

    %Class UID
    ClassUID = obj.rtPlanClassUID;
    meta.MediaStorageSOPClassUID = ClassUID;
    meta.SOPClassUID             = ClassUID;
    meta.FrameOfReferenceUID     = obj.FrameOfReferenceUID;
    %TransferSyntaxUID = '1.2.840.10008.1.2.1'; %Explicit VR Little Endian - correct?
    %meta.TransferSyntaxUID = TransferSyntaxUID;

    %Identifiers for this object
    try
        meta.SOPInstanceUID = obj.rtPlanMeta.SOPInstanceUID;   
    catch
        meta.SOPInstanceUID = dicomuid;
    end    
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
end

%Correct Positioning Offset
ct = obj.ct;
if ~any(isfield(ct,{'x','y','z'}))
    %positionOffset = transpose(ct.cubeDim ./ 2);
    %positionOffset = ct.cubeDim .* [ct.resolution.y, ct.resolution.x, ct.resolution.z] ./ 2;
    positionOffset = [ct.resolution.y, ct.resolution.x, ct.resolution.z] ./ 2;
    ct.x = ct.resolution.x*[0:ct.cubeDim(2)-1] - positionOffset(2);
    ct.y = ct.resolution.y*[0:ct.cubeDim(1)-1] - positionOffset(1);
    ct.z = ct.resolution.z*[0:ct.cubeDim(3)-1] - positionOffset(3);
end

positionOffsetCoordinates = [ct.x(1) ct.y(1) ct.z(1)];
positionOffsetCoordinates = positionOffsetCoordinates - [ct.resolution.x ct.resolution.y ct.resolution.z];



%Write references to image, RTStruct, RTDose
%RTStruct
meta.ReferencedStructureSetSequence.Item_1.ReferencedSOPClassUID = obj.rtStructClassUID;
meta.ReferencedStructureSetSequence.Item_1.ReferencedSOPInstanceUID = obj.rtssMeta.SOPInstanceUID;

%RTDose - %TODO the reference to the RTDose, also a little fishy because we
%need cross-references to the rtplan in rtdose to. Maybe we need to set-up
%IDs while constructing so we can set the references accordingly
%meta.DoseReferenceSquence.Item_1. ...


if obj.pln.propStf.numOfBeams ~= numel(obj.stf)
    matRad_cfg.error('Inconsistency in stf! number of beams not matching!');
end


%Sequences
%Sequences
%ToleranceTableSequence - Optional, we do not write this

% Fraction Sequence
%  meta.FractionGroupSequence.Item_1.
meta.FractionGroupSequence.Item_1.FractionGroupNumber =  1;
meta.FractionGroupSequence.Item_1.NumberOfFractionsPlanned =  obj.pln.numOfFractions;
meta.FractionGroupSequence.Item_1.NumberOfBeams = obj.pln.propStf.numOfBeams;
meta.FractionGroupSequence.Item_1.NumberOfBrachyApplicationSetups = 0;
% meta.FractionGroupSequence.Item_1.BeamDoseMeaning = 'FRACTION_LEVEL'; %TODO: This is probably no longer necessary.

refDoseSeq = struct();

%We need the doses to be exported already if we want to store the reference
for i = 1:numel(obj.rtDoseMetas)
    currItemStr = sprintf('Item_%d',i);
    refDoseSeq.(currItemStr).ReferencedSOPInstanceUID = obj.rtDoseMetas(i).SOPInstanceUID;
    refDoseSeq.(currItemStr).ReferencedSOPClassUID = obj.rtDoseMetas(i).SOPClassUID;
end

if ~isempty(refDoseSeq)
    meta.ReferencedDoseSequence = refDoseSeq;
    meta.FractionGroupSequence.Item_1.ReferencedDoseSequence = refDoseSeq;
end


for iBeam = 1:obj.pln.propStf.numOfBeams
    matRad_cfg.dispInfo('\tBeam %d: ',iBeam);
    
    %PatientSetupSequence - Required
    currBeamItemStr = sprintf('Item_%d',iBeam);
    meta.PatientSetupSequence.(currBeamItemStr) = struct('PatientPosition','HFS','PatientSetupNumber',iBeam,'SetupTechniqueDescription','');


    
    % Write Photon or ion RTPLAN
    if strcmp(obj.pln.radiationMode,'photons')
        BeamParam = 'BeamSequence';
        ControlParam = 'ControlPointSequence';
    elseif strcmp(obj.pln.radiationmode, 'protons') ||strcmp(obj.pln.radiationmode, 'helium') || strcmp(obj.pln.radiationmode, 'carbon')
        BeamParam = 'IonBeamSequence';
        ControlParam = 'IonControlPointSequence';
    else
        matRad_cfg.DispError('Not supported radiation mode of DICOM RT plan file.')
    end

    if strcmp(obj.pln.radiationMode,'photons')

        %seqItem references to the current BeamSequence
        beamSeqItem.TreatmentMachineName = obj.pln.machine;
        beamSeqItem.PrimaryDosimeterUnit = 'MU';
        beamSeqItem.SourceAxisDistance = obj.stf(iBeam).SAD;
        beamSeqItem.BeamNumber = iBeam;
        beamSeqItem.BeamType = 'DYNAMIC';
        beamSeqItem.RadiationType = 'PHOTON';
        beamSeqItem.TreatmentDeliveryType = 'TREATMENT';
        beamSeqItem.NumberOfWedges = 0;
        beamSeqItem.NumberOfCompensators = 0;
        beamSeqItem.NumberOfBoli = 0;
        beamSeqItem.NumberOfBlocks = 0;

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

        %Note also that if Final Cumulative Meterset Weight (300A,010E)
        %is equal to 100, then Cumulative Meterset Weight (300A,0134)
        %becomes equivalent to the percentage of Beam Meterset
        %(300A,0086) delivered at each control point. If Final
        %Cumulative Meterset Weight (300A,010E) is equal to Beam
        %Meterset (300A,0086), then the Cumulative Meterset Weight
        %(300A,0134) at each control point becomes equal to the
        %cumulative Meterset delivered at that control point.

        %This means that this can all be relative, when this value is
        %100, but we need an absolute value somewhere else. So let's
        %just use the shape weights and sum them up in the end
        %beamSeqItem.FinalCumulativeMetersetWeight = 100;


        
        % Take the information from the apertureInfo field
        %TODO: Get info about jaws? 
        if ~isfield(obj.resultGUI,'apertureInfo')
            matRad_cfg.dispError('Sequenced Apertures not found!')
        end
        apertureInfo = obj.resultGUI.apertureInfo;
        currBeamApertures = apertureInfo.beam(iBeam);

        %Get basic information about the collimator:
        limitDeviceSeq = struct();
        limitDeviceSeq.Item_1.RTBeamLimitingDeviceType = 'MLCX';        
        limitDeviceSeq.Item_1.NumberOfLeafJawPairs = apertureInfo.numOfMLCLeafPairs;

        if ~isfield(apertureInfo,'leafBoundaries') %TODO Information should be added by the sequencers
            matRad_cfg.dispWarning('Leaf Boundaries not specified in aperture info!');
            centralLeafPair = currBeamApertures.centralLeafPair;
            leafBoundaries = 0:apertureInfo.bixelWidth:(apertureInfo.numOfMLCLeafPairs)*apertureInfo.bixelWidth;
            centralLeafOffset = leafBoundaries(centralLeafPair) + apertureInfo.bixelWidth/2;
            leafBoundaries = leafBoundaries - centralLeafOffset;
        else
            leafBoundaries = apertureInfo.leafBoundaries;
        end
        limitDeviceSeq.Item_1.LeafPositionBoundaries = leafBoundaries;

        if ~isfield(apertureInfo,'SCD') %TODO implement that sequencers add this information (from machine file?)
            load(fullfile(matRad_cfg.matRadRoot,'basedata',['photons_' obj.pln.machine '.mat']),'machine');
            SCD = machine.meta.SCD;
        else
            SCD = aperturmateInfo.SCD;
        end
        limitDeviceSeq.Item_1.SourceToBeamLimitingDeviceDistance = SCD;
        
        beamSeqItem.BeamLimitingDeviceSequence = limitDeviceSeq;

        
        beamSeqItem.NumberOfControlPoints = currBeamApertures.numOfShapes;

        %Get IMRT shapes. Only IMRT support (DYNAMIC), no VMAT
        cumulativeMetersetWeight = 0;
        for iShape = 1:currBeamApertures.numOfShapes
            currCtrlSeqItemStr = sprintf('Item_%d',iShape);
            currCtrlSeqItem = struct();
            currCtrlSeqItem.ControlPointIndex = iShape-1;

            currShape = currBeamApertures.shape(iShape);

            %Write static Attributes only required for first control point:
            if iShape == 1
                currCtrlSeqItem.GantryAngle = obj.stf(iBeam).gantryAngle;
                currCtrlSeqItem.GantryRotationDirection = 'NONE';
                currCtrlSeqItem.PatientSupportAngle = obj.stf(iBeam).couchAngle;
                currCtrlSeqItem.PatientSupportRotationDirection = 'NONE';
                currCtrlSeqItem.TableTopEccentricAngle = 0;
                currCtrlSeqItem.TableTopEccentricRotationDirection = 'NONE';
                currCtrlSeqItem.IsocenterPosition = (obj.stf(iBeam).isoCenter + positionOffsetCoordinates)';
            end

            %Check energy, currently can only be constant
            energy = unique([obj.stf(iBeam).ray.energy]);
            if numel(energy) > 1
                matRad_cfg.dispError('Multiple Energies currently not supported!');
            end
            currCtrlSeqItem.NominalBeamEnergy = energy;

            %Collimator (& Jaws, but we leave them out for now)
            limitPosSeq = struct();
            limitPosSeq.Item_1.RTBeamLimitingDeviceType = 'MLCX';

            leftLeafPos = zeros(apertureInfo.numOfMLCLeafPairs,1);
            rightLeafPos = leftLeafPos;

            leftLeafPos(logical(currBeamApertures.isActiveLeafPair)) = currShape.leftLeafPos;
            rightLeafPos(logical(currBeamApertures.isActiveLeafPair)) = currShape.rightLeafPos;

            limitPosSeq.Item_1.LeafJawPositions = [leftLeafPos; rightLeafPos];

            currCtrlSeqItem.BeamLimitingDevicePositionSequence = limitPosSeq;

            %This always needs to be zero for the first Item. However,
            %the last one should be the fully cumulated value. Does
            %this mean that the last control point ?
            currCtrlSeqItem.CumulativeMetersetWeight = cumulativeMetersetWeight;

            shapeWeight = currBeamApertures.shape(iShape).weight; %Is this correct? It basically means that the next
            cumulativeMetersetWeight = cumulativeMetersetWeight + shapeWeight;

            beamSeqItem.ControlPointSequence.(currCtrlSeqItemStr) = currCtrlSeqItem;
        end

        %We need to add a final control point only containing the
        %cumulative meterset
        currCtrlSeqItemStr = sprintf('Item_%d',iShape+1);
        beamSeqItem.ControlPointSequence.(currCtrlSeqItemStr) = struct();
        beamSeqItem.ControlPointSequence.(currCtrlSeqItemStr).CumulativeMetersetWeight = cumulativeMetersetWeight;
        
        %Add the final cumulative meterset weight (see comment above, we
        %use absolute values and not 100, so we can easy add it to the Fraction sequence)
        beamSeqItem.FinalCumulativeMetersetWeight = cumulativeMetersetWeight;

        %TODO: Do we need anything else here?   

        matRad_cfg.dispInfo('exported %d Control Points.\n',iShape+1);
    end

    meta.BeamSequence.(currBeamItemStr) = beamSeqItem;

    %TODO meta.FractionGroupSequence.Item_1.ReferencedBeamSequence
    refBeamSeqItem = struct();
    resultGUIbeamStr = sprintf('physicalDose_beam%d',iBeam);
    if isfield(obj.resultGUI,resultGUIbeamStr)
        refBeamSeqItem.BeamDose = mean(obj.resultGUI.physicalDose_beam1,"all"); %TODO: Probably no longer necessary according to standard
    end
    refBeamSeqItem.BeamMeterset = cumulativeMetersetWeight;
    refBeamSeqItem.ReferencedBeamNumber = iBeam;
    meta.FractionGroupSequence.Item_1.ReferencedBeamSequence.(currBeamItemStr) = refBeamSeqItem;
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
obj.rtPlanMeta = meta;

end
