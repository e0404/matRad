function [cube, metadata] = matRad_writeNifTI(filepath,cube,metadata)
% matRad NifTI reader
% 
% call
%   [cube, metadata] = matRad_writeNifTI(filename)
%
% input
%   filename:   full path to .nii(.gz) file
%
% output
%   file will be written to disk
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%Instantiate matrad config
matRad_cfg = MatRad_Config.instance();

%Check for Image Processing toolbox
avail = license('test', 'image_toolbox');
if ~avail
    matRad_cfg.dispError('Image Processing Toolbox is required for reading NifTI files!\n');
else
    success = license('checkout', 'image_toolbox');
    if ~success
        matRad_cfg.dispError('Image Processing Toolbox license could not be checked out!\n');
    end
end

info = struct();

[~,~,endian] = computer;
switch endian
    case 'L'
        endian = 'little';
    case 'B'
        endian = 'big';
    otherwise
        matRad_cfg.dispError('Unknown endian!');
end

%Check if compression was defined
if ~isfield(metadata,'compress') || ~islogical(metadata.compress)
    metadata.compress = true;
end  

%Default Matlab Axis Permuation
metadata.axisPermutation = [2 1 3];
cube = permute(cube,metadata.axisPermutation);

info.PixelDimensions = [metadata.resolution(metadata.axisPermutation)];
info.Datatype = metadata.datatype;
info.ImageSize = size(cube);
info.Description = sprintf('Exported from matRad %s',matRad_version());
if length(info.Description) >= 80
    info.Description = info.Description(1:79);
end

if max(info.ImageSize) > 32767
    info.Version = 'NIfTI2';
else
    info.Version = 'NIfTI1';
end

info.Qfactor = 1;
info.SpaceUnits = 'Millimeter';
info.TimeUnits = 'None';
info.SliceCode = 'Unknown';
info.FrequencyDimension = 0;
info.PhaseDimension = 0;
info.SpatialDimension = 0;

%Set up Transform Matrix
info.TransformName = 'Sform';

T=eye(4);


%Correct for coordinate system
switch metadata.coordinateSystem
    case 'LPS'
        tmpT = eye(4,4);
        tmpT(1:2,1:2) = -1*tmpT(1:2,1:2);
        T = T*tmpT;
    otherwise
        matRad_cfg.dispError('Only LPS currently supported for export!');
end

%Add resolution
T(:,1) = T(:,1) * metadata.resolution(metadata.axisPermutation(1));
T(:,2) = T(:,2) * metadata.resolution(metadata.axisPermutation(2));
T(:,3) = T(:,3) * metadata.resolution(metadata.axisPermutation(3));

%Now add Translation
T(4,1) = metadata.imageOrigin(metadata.axisPermutation(1));
T(4,2) = metadata.imageOrigin(metadata.axisPermutation(2));
T(4,3) = metadata.imageOrigin(metadata.axisPermutation(3));

info.Transform = affine3d(T);
cube = cast(cube,metadata.datatype);

niftiwrite(cube,filepath,info,'Endian',endian,'Compressed',metadata.compress);



