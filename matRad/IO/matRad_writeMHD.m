function matRad_writeMHD(filepath,cube,metadata)
% matRad function to write mha files
% 
% call
%   matRad_writeMHA(filepath,cube,metadata)
%
% input
%   filepath:   full filename (with extension)
%   cube:       3D array to be written into file
%   metadata:   struct of metadata. Writer will wrap the existing metadata 
%               to MHA standard-specific fields 
%               Necessary fieldnames are:
%               - resolution: [x y z]
%               - datatype: numeric MATLAB-Datatype
%
% output
%   file will be written to disk
%
% References
%   https://itk.org/Wiki/MetaIO/Documentation
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

%Sanity checks and restrictions
dimensions = size(cube);
if numel(dimensions) ~= 3
    matRad_cfg.dispError('Sorry! matRad only supports 3-dimensional MHD output');
end

fid = fopen(filepath, 'wb');
if fid <= 0
    matRad_cfg.dispError('Could not open MHD destination file!');
end

%We perform the permutation
if ~isfield(metadata, 'axisPermutation')
    % This reverts the matRlab conventianl indexing
    axisPermutation = [2,1,3];
else
    if ~isequal(metadata.axisPermutation, [2,1,3])
        matRad_cfg.dispWarning('Unconventianal permutation of patient indexing, this might cause inconsistency');
    end
    axisPermutation = metadata.axisPermutation;
end

% Force the permutation here according to the axis permutation
cube = permute(cube, axisPermutation);

% Need to permute the dimensions as well
dimensions = size(cube);

% Note in the cube permutation:
% Permutation of the cube is enforced here to standdard indexing so that
% reconstruction of the cube for further use does not rely on the use of 
% a transformation matrix, which is now only the identity matrix.

%The transformation matrix is now the unit matrix
transformMatrix = diag(ones(1,numel(dimensions)));
tmString = sprintf(' %d',transformMatrix(:));

%Determine the endian
[~,~,endian] = computer;
switch endian
    case 'L'
        byteOrderMSB = 'False';
    case 'B'
        byteOrderMSB = 'True';
    otherwise
        error('Unknown endian!');
end

[path,name,ext] = fileparts(filepath);
filenameRaw = [name '.raw'];

fprintf(fid, 'ObjectType = Image\n');
fprintf(fid, 'NDims = %d\n',numel(dimensions));
fprintf(fid, 'BinaryData = True\n');
fprintf(fid, 'BinaryDataByteOrderMSB = %s\n',byteOrderMSB); %Not sure about this field
%fprintf(fid, 'ElementByteOrderMSB = %s\n',byteOrderMSB); %Not sure about this field
fprintf(fid, 'TransformMatrix =%s\n',tmString);
fprintf(fid, 'Offset = %f %f %f\n',metadata.imageOrigin(1),metadata.imageOrigin(2),metadata.imageOrigin(3));
fprintf(fid, 'AnatomicalOrientation = RAI\n'); %Did not double check this line
fprintf(fid, 'ElementSpacing = %f %f %f\n',metadata.resolution(1),metadata.resolution(2),metadata.resolution(3));
fprintf(fid, 'DimSize = %d %d %d\n',dimensions(1),dimensions(2),dimensions(3));
fprintf(fid, 'ElementType = %s\n',matlabTypeToMHAtype(metadata.datatype));
fprintf(fid, sprintf('ElementDataFile = %s\n',filenameRaw));
fclose(fid);

%% write data file
dataFileHandle = fopen(fullfile(path,filenameRaw),'w');
fwrite(dataFileHandle,cube(:),metadata.datatype,lower(endian));
fclose(dataFileHandle);

end

function newType = matlabTypeToMHAtype(datatype)
switch datatype
    case {'single','float'}
        newType = 'MET_FLOAT';     
    case 'double'
        newType = 'MET_DOUBLE';
    case {'uchar','uint8'}
       newType = 'MET_UCHAR';
    case {'logical','int8','char'}
        newType = 'MET_CHAR';
    case 'int16'
        newType = 'MET_SHORT';
    case 'uint16'
        newType = 'MET_USHORT';
    case 'int32'
        newType = 'MET_INT';
    case 'int64'
        newType = 'MET_LONG';
    case 'uint32'
        newType = 'MET_UINT';
    case 'uint64'
        newType = 'MET_ULONG';
    otherwise
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError(['Datatype ' datatype ' not supported by MHD exporter!']);
end
end



