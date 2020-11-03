function pln = matRad_importDicomRTPlan(ct, rtPlanFiles, dicomMetaBool)
% matRad function to import dicom RTPLAN data
% 
% call
%   pln = matRad_importDicomRTPlan(ct, rtPlanFiles, dicomMetaBool)
%
% input
%   ct:             ct imported by the matRad_importDicomCt function
%   rtPlanFiles:   	list of RTPlan Dicom files
%   dicomMetaBool:  import whole dicom information
%
% output
%   pln:            matRad pln struct with meta information. Note that
%                   bixelWidth is determined via the importSteering function.
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

%% load plan file
% check size of RT Plan
if size(rtPlanFiles,1) ~= 1
   errordlg('Too few or to many RTPlan files')
end

% read information out of the RT file
if verLessThan('matlab','9')
    planInfo = dicominfo(rtPlanFiles{1});
else
    planInfo = dicominfo(rtPlanFiles{1},'UseDictionaryVR',true);
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
end

% transform iso. At the moment just this way for HFS
if ct.dicomInfo.ImageOrientationPatient == [1;0;0;0;1;0]
    isoCenter = isoCenter - ones(length(BeamSeqNames),1) * ...
        ([ct.x(1) ct.y(1) ct.z(1)] - [ct.resolution.x ct.resolution.y ct.resolution.z]);
else
    error('This Orientation is not yet supported.');
end

%% read constant parameters
% readout charge and mass to set radiationMode to matRad specific name
radiationMode = planInfo.(BeamParam).Item_1.RadiationType;
if ~strncmpi(radiationMode,'photons',6)
    try
        radiationMass = planInfo.(BeamParam).Item_1.RadiationMassNumber;
        radiationAtomicNumber = planInfo.(BeamParam).Item_1.RadiationAtomicNumber;
    catch
        warning('Could not determine mass and atomic number of the particle');
    end
end

if strncmpi(radiationMode,'photons',6)
    radiationMode = 'photons';
elseif strncmpi(radiationMode,'proton',6)
    radiationMode = 'protons';
elseif (strncmpi(radiationMode,'ion',3) && radiationMass == 12 && radiationAtomicNumber == 6)
    radiationMode = 'carbon';
else
    warning('The given type of radiation is not yet supported');
end

% extract field shapes
if strcmp(radiationMode, 'photons')
           
    fractionSequence         = planInfo.FractionGroupSequence.Item_1;
    pln.propStf.collimation  = matRad_importFieldShapes(BeamSequence,fractionSequence);
    
end

%% write parameters found to pln variable
pln.radiationMode   = radiationMode; % either photons / protons / carbon
pln.numOfFractions  = planInfo.FractionGroupSequence.Item_1.NumberOfFractionsPlanned;
pln.machine         = planInfo.(BeamParam).Item_1.TreatmentMachineName;

pln.propStf.isoCenter    = isoCenter;
pln.propStf.bixelWidth   = NaN; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles = [gantryAngles{1:length(BeamSeqNames)}];
pln.propStf.couchAngles  = [PatientSupportAngle{1:length(BeamSeqNames)}]; % [°]
pln.propStf.numOfBeams   = length(BeamSeqNames);


pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.propOpt.runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles

% if we imported field shapes then let's trigger field based dose calc by
% setting the bixelWidth to 'field'
if isfield(pln.propStf,'collimation')
    pln.propStf.bixelWidth  = 'field'; 
end

% timestamp
pln.DicomInfo.timeStamp = datestr(clock);

try
   pln.DicomInfo.SOPClassUID = planInfo.SOPClassUID;
   pln.DicomInfo.SOPInstanceUID = planInfo.SOPInstanceUID;
   pln.DicomInfo.ReferencedDoseSequence = planInfo.ReferencedDoseSequence;
catch
end

% safe entire dicomInfo
if dicomMetaBool == true
    pln.DicomInfo.Meta = planInfo;
end
end
