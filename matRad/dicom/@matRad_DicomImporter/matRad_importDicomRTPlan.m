function obj = matRad_importDicomRTPlan(obj)
% matRad function to import dicom RTPLAN data
% 
% In your object, there must be properties that contain:
%   - ct imported by the matRad_importDicomCt function;
%   - list of RTPlan Dicom files;
%   - a boolean, if you don't want to import whole dicom information set it
%   false.
% 
% Output - matRad pln structure with meta information.
% Note that bixelWidth is determined via the importSteering function.
%
% call
%   obj = matRad_importDicomRTPlan(obj)
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
matRad_checkEnvDicomRequirements(matRad_cfg.env);

%% load plan file
% check size of RT Plan
if size(obj.importFiles.rtplan,1) ~= 1
   errordlg('Too few or to many RTPlan files')
end

% read information out of the RT file
if matRad_cfg.isOctave || verLessThan('matlab','9')
    planInfo = dicominfo(obj.importFiles.rtplan{1});
else
    planInfo = dicominfo(obj.importFiles.rtplan{1},'UseDictionaryVR',true);
end

% check which type of Radiation is used
if isfield(planInfo, 'BeamSequence')
    BeamParam = 'BeamSequence';
    ControlParam = 'ControlPointSequence';
elseif isfield(planInfo, 'IonBeamSequence')
    BeamParam = 'IonBeamSequence';
    ControlParam = 'IonControlPointSequence';
else
    errordlg('Not supported kind of DICOM RT plan file.');    
end

% get beam sequence
BeamSequence = planInfo.(BeamParam);
BeamSeqNames = fieldnames(BeamSequence);

% use the treatment beams only
if isfield(BeamSequence.(BeamSeqNames{1}),'TreatmentDeliveryType')
    for i = 1:length(BeamSeqNames)
        currBeamSeq = BeamSequence.(BeamSeqNames{i});
        try
            treatDelType = currBeamSeq.TreatmentDeliveryType;
            if ~strcmpi(treatDelType,'TREATMENT')
                BeamSequence = rmfield(BeamSequence,BeamSeqNames{i});
            end
        catch
            warning('Something went wrong while determining the type of the beam.');
        end
    end
    BeamSeqNames = fieldnames(BeamSequence);
end


%% get information may change between beams
% loop over beams
gantryAngles{length(BeamSeqNames)} = [];
PatientSupportAngle{length(BeamSeqNames)} = [];
isoCenter = NaN*ones(length(BeamSeqNames),3);
for i = 1:length(BeamSeqNames)   
    currBeamSeq             = BeamSequence.(BeamSeqNames{i});
    % parameters not changing are stored in the first ControlPointSequence
    gantryAngles{i}         = currBeamSeq.(ControlParam).Item_1.GantryAngle;
    PatientSupportAngle{i}  = currBeamSeq.(ControlParam).Item_1.PatientSupportAngle;
        isoCenter(i,:)          = currBeamSeq.(ControlParam).Item_1.IsocenterPosition';
    if isoCenter(i,1) < min(obj.ct.x) || isoCenter(i,1) > max(obj.ct.x)  || isoCenter(i,2) < min(obj.ct.y) || ...
            isoCenter(i,2) > max(obj.ct.y) || isoCenter(i,3) < min(obj.ct.z) || isoCenter(i,3) > max(obj.ct.z)
        isoCenter(i,:)          = matRad_getIsoCenter(obj.cst, obj.ct);
    end
end

% transform iso. At the moment just this way for HFS
if obj.ct.dicomInfo.ImageOrientationPatient ~= [1;0;0;0;1;0]    
    matRad_cfg.dispError('This Orientation is not yet supported.');
end

%% read constant parameters
% readout charge and mass to set radiationMode to matRad specific name
radiationMode = planInfo.(BeamParam).Item_1.RadiationType;
if ~strncmpi(radiationMode,'photons',6)
    try
        radiationMass = planInfo.(BeamParam).Item_1.RadiationMassNumber;
        radiationAtomicNumber = planInfo.(BeamParam).Item_1.RadiationAtomicNumber;
    catch
        matRad_cfg.dispWarning('Could not determine mass and atomic number of the particle');
    end
end

if strncmpi(radiationMode,'photons',6)
    radiationMode = 'photons';
elseif strncmpi(radiationMode,'proton',6)
    radiationMode = 'protons';
elseif (strncmpi(radiationMode,'ion',3) && radiationMass == 12 && radiationAtomicNumber == 6)
    radiationMode = 'carbon';
else
    matRad_cfg.dispError('The given type of radiation is not yet supported');
end

% extract field shapes
if strcmp(radiationMode, 'photons')
           
    fractionSequence         = planInfo.FractionGroupSequence.Item_1;
    obj.pln.propStf.collimation  = matRad_importFieldShapes(BeamSequence,fractionSequence);
    
end

%% write parameters found to pln variable
obj.pln.radiationMode   = radiationMode; % either photons / protons / carbon
obj.pln.numOfFractions  = planInfo.FractionGroupSequence.Item_1.NumberOfFractionsPlanned;

% set handling of multiple scenarios -> default: only nominal
obj.pln.multScen = matRad_multScen(obj.ct,'nomScen');
if isfield(BeamSequence.Item_1, 'TreatmentMachineName')
    obj.pln.machine         = BeamSequence.Item_1.TreatmentMachineName;
else 
    obj.pln.machine         = 'Generic';
end
% set bio model parameters (default physical opt, no bio model)
obj.pln.bioModel = matRad_bioModel(obj.pln.radiationMode,'none');

% set properties for steering
obj.pln.propStf.isoCenter    = isoCenter;
obj.pln.propStf.bixelWidth   = NaN; % [mm] / also corresponds to lateral spot spacing for particles
obj.pln.propStf.gantryAngles = [gantryAngles{1:length(BeamSeqNames)}];
obj.pln.propStf.couchAngles  = [PatientSupportAngle{1:length(BeamSeqNames)}]; % [??]
obj.pln.propStf.numOfBeams   = length(BeamSeqNames);
numOfVoxels = 1;
for i = 1:length(obj.ct.cubeDim)
    numOfVoxels = numOfVoxels*obj.ct.cubeDim(i);
end
obj.pln.numOfVoxels          = numOfVoxels;
obj.pln.VoxelDimentions      = obj.ct.cubeDim;

%if there is not special doseGrid for rtdose
if ~obj.importFiles.useImportGrid && isfield(obj.importFiles,'rtdose')
    obj.pln.propDoseCalc.doseGrid.resolution.x = obj.ct.resolution.x;
    obj.pln.propDoseCalc.doseGrid.resolution.y = obj.ct.resolution.y;
    obj.pln.propDoseCalc.doseGrid.resolution.z = obj.ct.resolution.z;
end

% turn off sequerncing an DAO by default
obj.pln.propOpt.runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
obj.pln.propOpt.runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles

% if we imported field shapes then let's trigger field based dose calc by
% setting the bixelWidth to 'field'
if isfield(obj.pln.propStf,'collimation')
    obj.pln.propStf.bixelWidth  = 'field'; 
end

% timestamp
obj.pln.DicomInfo.timeStamp = datestr(clock);

try
   obj.pln.DicomInfo.SOPClassUID = planInfo.SOPClassUID;
   obj.pln.DicomInfo.SOPInstanceUID = planInfo.SOPInstanceUID;
   obj.pln.DicomInfo.ReferencedDoseSequence = planInfo.ReferencedDoseSequence;
catch
end

% safe entire dicomInfo
if obj.dicomMetaBool == true
    obj.pln.DicomInfo.Meta = planInfo;
end
end
