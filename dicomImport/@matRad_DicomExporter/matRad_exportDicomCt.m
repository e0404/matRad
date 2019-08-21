function obj = matRad_exportDicomCt(obj)
% matRad function to export dicom ct data
%
% call
%   ct = matRad_exportDicomCt(ctList, resolution, dicomMetaBool, visBool)
%
% input
%   dicomDir:       Directory to store dicom files
%   ct:             matRad ct struct
%   meta:           pre-given dicom meta data
%
% output
%   status:         status from the dicomwrite function as struct array
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

%default meta
meta.PatientID = obj.PatientID;
meta.PatientName = obj.PatientName;
meta.PatientPosition = obj.PatientPosition;
meta.StudyID = obj.StudyID;
meta.StudyDate = obj.StudyDate;
meta.StudyTime = obj.StudyTime;
meta.StudyInstanceUID = obj.StudyInstanceUID;
meta.FrameOfReferenceUID = obj.FrameOfReferenceUID;

obj.ctMeta.SeriesInstanceUID = dicomuid;
meta.SeriesInstanceUID = obj.ctMeta.SeriesInstanceUID;

ct = obj.ct;

nSlices = ct.cubeDim(3);
%Create X Y Z vectors if not present
if ~any(isfield(ct,{'x','y','z'}))
    positionOffset = transpose(ct.cubeDim ./ 2);
    ct.x = ct.resolution.x*[0:ct.cubeDim(1)-1] - positionOffset;
    ct.y = ct.resolution.y*[0:ct.cubeDim(2)-1] - positionOffset;
    ct.z = ct.resolution.z*[0:ct.cubeDim(3)-1] - positionOffset;
end

obj.ct = ct;

%Since we are exporting HU directly --> no rescaling in any case
meta.SliceThickness = ct.resolution.z;
meta.PixelSpacing = [ct.resolution.y; ct.resolution.x];
meta.ImageOrientationPatient = [1;0;0;0;1;0]; %lps
meta.RescaleType = 'HU';

if isfield(ct,'z')
    z = ct.z;
end

ctCube = ct.cubeHU{1};
ctMin = min(ctCube(:));
ctCube = ctCube - ctMin;
ctMax = max(ctCube(:));
ctCube = ctCube ./ ctMax;
ctCube = uint16(ctCube*double(uint16(Inf)));

meta.RescaleIntercept = ctMin;
meta.RescaleSlope = ctMax / double(uint16(Inf) + 1);

meta.ImageType = 'DERIVED\PRIMARY\AXIAL';

fileName = 'ct_slice_';

obj.ctSliceMetas = struct([]);
obj.ctExportStatus = struct([]);

for i = 1:nSlices
    ctSlice = ctCube(:,:,i);
    %ctSlice = permute(ctSlice,[1 2]);
    
    obj.ctSliceMetas = obj.addStruct2StructArray(obj.ctSliceMetas,meta);
    
    obj.ctSliceMetas(i).ImagePositionPatient = [ct.x(1); ct.y(1); ct.z(i)];
    obj.ctSliceMetas(i).SlicePosition = z(i);
    
    %Create and store unique ID
    obj.ctSliceMetas(i).SOPInstanceUID = dicomuid;
    obj.ctSliceMetas(i).SOPClassUID = '1.2.840.10008.5.1.4.1.1.2';
    obj.ctSliceMetas(i).MediaStorageSOPInstanceUID = obj.ctSliceMetas(i).SOPInstanceUID;
    
    fullFileName = fullfile(obj.dicomDir,[fileName num2str(i) '.dcm']);
    
     status = dicomwrite(ctSlice,fullFileName,obj.ctSliceMetas(i),'ObjectType','CT Image Storage');
     obj.ctExportStatus = obj.addStruct2StructArray(obj.ctExportStatus,status);
end

end

