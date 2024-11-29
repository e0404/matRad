function obj = matRad_interpDicomDoseCube(obj) 
% matRad function to interpolate a given Dicom Dose Cube dicom RTDOSE data
%
% In your object, there must be properties that contain:
%   - ct imported by the matRad_importDicomCt function;
%   - one (of several) dose cubes which should be interpolated.
% Optional:
%   - pln structure.
%
% Output - structure with different actual current dose cube and several 
% meta data.
%
% call
%   obj = matRad_interpDicomDoseCube(obj)
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

% read information out of the RT file

matRad_cfg = MatRad_Config.instance();
matRad_checkEnvDicomRequirements(matRad_cfg.env);

dosefile = obj.importRTDose.currDose{1};

if matRad_cfg.isOctave || verLessThan('matlab','9')
    doseInfo = dicominfo(dosefile);
else
    doseInfo = dicominfo(dosefile,'UseDictionaryVR',true);
end

% read the dosefile itself
dosedata = dicomread(dosefile);
obj.importRTDose.dose.cube = double(dosedata);

% give it an internal name
 % obj.dose.internalName = obj.currDose{12};%?????

% read out the resolution
obj.importRTDose.dose.resolution.x = doseInfo.PixelSpacing(1);
obj.importRTDose.dose.resolution.y = doseInfo.PixelSpacing(2);
obj.importRTDose.dose.resolution.z = obj.importFiles.resz;

% target resolution is ct.resolution
target_resolution = obj.ct.resolution;

% convert dosedata to 3-D cube
obj.importRTDose.dose.cube = squeeze(obj.importRTDose.dose.cube(:,:,1,:));

% ct resolution is target resolution, now convert to new cube;

% generating grid vectors
x = doseInfo.ImagePositionPatient(1) + doseInfo.ImageOrientationPatient(1) * ...
                                       doseInfo.PixelSpacing(1) * double(0:doseInfo.Columns - 1);
y = doseInfo.ImagePositionPatient(2) + doseInfo.ImageOrientationPatient(5) * ...
                                       doseInfo.PixelSpacing(2) * double(0:doseInfo.Rows - 1);
z = doseInfo.ImagePositionPatient(3) + doseInfo.GridFrameOffsetVector;

% set up grid matrices - implicit dimension permuation (X Y Z-> Y X Z)
% Matlab represents internally in the first matrix dimension the
% ordinate axis and in the second matrix dimension the abscissas axis
[ X,  Y,  Z] = meshgrid(x,y,z);
[Xq, Yq, Zq] = meshgrid(obj.ct.x,obj.ct.y,obj.ct.z);

% get GridScalingFactor
gridScale = double(doseInfo.DoseGridScaling);
% rescale importRTDose.dose.cube
obj.importRTDose.dose.cube = gridScale * obj.importRTDose.dose.cube;

% interpolation to ct grid - cube is now stored in Y X Z
obj.importRTDose.dose.cube = interp3(X,Y,Z,obj.importRTDose.dose.cube,Xq,Yq,Zq,'linear',0);

% write new parameters
obj.importRTDose.dose.resolution = obj.ct.resolution;
obj.importRTDose.dose.x = obj.ct.x;
obj.importRTDose.dose.y = obj.ct.y;
obj.importRTDose.dose.z = obj.ct.z;

% write Dicom-Tags
obj.importRTDose.dose.dicomInfo.PixelSpacing            = [target_resolution.x; ...
                                                target_resolution.y];
obj.importRTDose.dose.dicomInfo.ImagePositionPatient    = [min(obj.importRTDose.dose.x); min(obj.importRTDose.dose.y); min(obj.importRTDose.dose.z)];
obj.importRTDose.dose.dicomInfo.SliceThickness          = target_resolution.z;
obj.importRTDose.dose.dicomInfo.ImageOrientationPatient = doseInfo.ImageOrientationPatient;
obj.importRTDose.dose.dicomInfo.DoseType                = doseInfo.DoseType;
obj.importRTDose.dose.dicomInfo.DoseSummationType       = doseInfo.DoseSummationType;
%importRTDose.dose.dicomInfo.InstanceNumber          = doseInfo.InstanceNumber; %Not
%always given
obj.importRTDose.dose.dicomInfo.SOPClassUID             = doseInfo.SOPClassUID;
obj.importRTDose.dose.dicomInfo.SOPInstanceUID          = doseInfo.SOPInstanceUID;
if isfield(doseInfo,'ReferencedRTPlanSequence')
    obj.importRTDose.dose.dicomInfo.ReferencedRTPlanSequence = doseInfo.ReferencedRTPlanSequence;
end

end

