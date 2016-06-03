function pln = matRad_importDicomRTPlan(ct, rtPlanFiles, dicomMetaBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to import dicom RTPLAN data
% 
% call
%   pln = matRad_importDicomRTPlan(ct, rtPlanFiles, dicomMetaBool)
%
% input
%   ct:             ct imported by the matRad_importDicomCt function
%   rtDoseFiles:   	list of RTDOSE Dicom files
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


% %% print current status of the import script
% fprintf('Please provide an appropriate machine file:\n');

%% load plan file
% check size of RT Plan
if size(rtPlanFiles,1) ~= 1
   errordlg('Too few or to many RTPlan files')
end

% read information out of the RT file
planInfo = dicominfo(rtPlanFiles{1});

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

% use the treatment beams only
BeamSequence = planInfo.(BeamParam);
BeamSeqNames = fieldnames(BeamSequence);
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

%% get information may change between beams
% loop over beams
gantryAngles{length(BeamSeqNames)} = [];
PatientSupportAngle{length(BeamSeqNames)} = [];
isoCenter{length(BeamSeqNames)} = [];
for i = 1:length(BeamSeqNames)   
    currBeamSeq             = BeamSequence.(BeamSeqNames{i});
    % parameters not changing are stored in the first ControlPointSequence
    gantryAngles{i}         = currBeamSeq.(ControlParam).Item_1.GantryAngle;
    PatientSupportAngle{i}  = currBeamSeq.(ControlParam).Item_1.PatientSupportAngle;
    isoCenter{i}            = currBeamSeq.(ControlParam).Item_1.IsocenterPosition;
end

% check wether isocenters are consistent
if numel(isoCenter) > 1
    if ~isequal(isoCenter{:})
       errordlg('Values for isocenter are not consistent.')
    end
end

% transform iso. At the moment just this way for HFS
if ct.dicomInfo.ImageOrientationPatient == [1;0;0;0;1;0]
    isoCenter = isoCenter{1}' - ct.dicomInfo.ImagePositionPatient' + ...
                         [ct.resolution.x ct.resolution.y ct.resolution.z];
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

%% write parameters found to pln variable
pln.isoCenter       = isoCenter;
pln.radiationMode   = radiationMode; % either photons / protons / carbon
pln.bixelWidth      = NaN; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [gantryAngles{1:length(BeamSeqNames)}];
pln.couchAngles     = [PatientSupportAngle{1:length(BeamSeqNames)}]; % [°]
pln.numOfBeams      = length(BeamSeqNames);
pln.numOfVoxels     = numel(ct.cube{1});
pln.voxelDimensions = ct.cubeDim;
pln.bioOptimization = NaN; % none: physical optimization; effect: effect-based optimization; RBExD: optimization of RBE-weighted dose
pln.numOfFractions  = planInfo.FractionGroupSequence.Item_1.NumberOfFractionsPlanned;
pln.runSequencing   = NaN; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = NaN; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.machine         = 'unknown';

% timestamp
pln.timeStamp = datestr(clock);

try
   pln.SOPClassUID = planInfo.SOPClassUID;
   pln.SOPInstanceUID = planInfo.SOPInstanceUID;
   pln.ReferencedDoseSequence = planInfo.ReferencedDoseSequence;
catch
end

% safe entire dicomInfo
if dicomMetaBool == true
    pln.dicomMeta = planInfo;
end
end
