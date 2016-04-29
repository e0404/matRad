function [ dose ] = matRad_interpDoseCube( ct, currDose )

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

% read information out of the RT file
dosefile = currDose{1};
info = dicominfo(dosefile);

% read the dosefile itself
dosedata = dicomread(dosefile);
% convert 16-bit integer to double precision, therefore to 1 normalized
dosedata = im2double(dosedata);

% give it an internal name
dose.internalName = currDose{12};

% read out the resolution
dose.resolution.x = info.PixelSpacing(1);
dose.resolution.y = info.PixelSpacing(2);
dose.resolution.z = info.SliceThickness;

% target resolution is ct.resolution
target_resolution = ct.resolution;

% convert dosedata to 3-D cube
dose.cube = squeeze(dosedata(:,:,1,:));

% ct resolution is target resolution, now convert to new cube;

% generating grid vectors
x = info.ImagePositionPatient(1) + info.PixelSpacing(1) * double([0:info.Columns - 1]);
y = info.ImagePositionPatient(2) + info.PixelSpacing(2) * double([0:info.Rows - 1]');
z = [info.ImagePositionPatient(3) + info.GridFrameOffsetVector]';

% new vectors
xq = [min(ct.x) : target_resolution.x : max(ct.x)]';
yq = [min(ct.y) : target_resolution.y : max(ct.y)]';
zq =  min(ct.z) : target_resolution.z : max(ct.z);

% scale cube from relative (normalized) to absolute values
% need BitDepth
bitDepth = double(info.BitDepth);
% get GridScalingFactor
gridScale = double(info.DoseGridScaling);
% CAUTION: Only valid if data is converted via im2double
doseScale = (2 ^ bitDepth - 1) * gridScale;
% rescale dose.cube
dose.cube = doseScale * dose.cube;

% interpolation to ct grid
dose.cube = interp3(x,y,z,dose.cube,xq,yq,zq,'linear',0);

% write new parameters
dose.resolution = ct.resolution;
dose.x = xq';
dose.y = yq';
dose.z = zq;

% check whether grid position are the same as the CT grid positions are
if ~(isequal(dose.x,ct.x) && isequal(dose.y,ct.y) && isequal(dose.z,ct.z))
    errordlg('CT-Grid and Dose-Grid are still not the same');
end

% write Dicom-Tags
dose.dicomInfo.PixelSpacing            = [target_resolution.x; ... 
                                                target_resolution.y];
dose.dicomInfo.ImagePositionPatient    = [min(dose.x); min(dose.y); min(dose.z)];
% only if z is evenly spaced!
dose.dicomInfo.SliceThickness          = target_resolution.z;
dose.dicomInfo.ImageOrientationPatient = info.ImageOrientationPatient;
dose.dicomInfo.DoseType                = info.DoseType;
dose.dicomInfo.DoseSummationType       = info.DoseSummationType;
dose.dicomInfo.InstanceNumber          = info.InstanceNumber;
dose.dicomInfo.SOPClassUID             = info.SOPClassUID;
dose.dicomInfo.SOPInstanceUID          = info.SOPInstanceUID;
dose.dicomInfo.ReferencedRTPlanSequence = info.ReferencedRTPlanSequence;

end

