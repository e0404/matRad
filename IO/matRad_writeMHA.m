function matRad_writeMHA(filepath,cube,metadata)
% matRad function to write mha files
% 
% call
%   matRad_writeMHA(cube,metadata,filename)
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
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sanity checks and restrictions
dimensions = size(cube);
if numel(dimensions) ~= 3
    error('Sorry! matRad only supports 3-dimensional MHA output');
end

fid = fopen(filepath, 'wb');
if fid <= 0
    error('Could not open MHA destination file!');
end
cleaner = onCleanup(@() fclose(fid));

%We perform the permutation
if isfield(metadata,'axisPermutation')
    cube = permute(cube,metadata.axisPermutation);
end
%The transformation matrix is now the unit matrix
transformMatrix = diag(ones(1,numel(dimensions)));
tmString = sprintf(' %d',transformMatrix(:));

%Determine the endian
[~,~,endian] = computer;
switch endian
    case 'L'
        byteOrderMSB = 'True';
    case 'B'
        byteOrderMSB = 'False';
    otherwise
        error('Unknown endian!');
end

fprintf(fid, 'ObjectType = Image\n');
fprintf(fid, 'NDims = %d\n',numel(dimensions));
fprintf(fid, 'BinaryData = True\n');
fprintf(fid, 'BinaryDataByteOrderMSB = %s\n',byteOrderMSB); %Not sure about this field
fprintf(fid, 'ElementByteOrderMSB = %s\n',byteOrderMSB); %Not sure about this field
fprintf(fid, 'TransformMatrix =%s\n',tmString);
fprintf(fid, 'Offset = %f %f %f\n',metadata.imageOrigin(1),metadata.imageOrigin(2),metadata.imageOrigin(3));
fprintf(fid, 'AnatomicalOrientation = RAI\n'); %Did not double check this line
fprintf(fid, 'ElementSpacing = %f %f %f\n',metadata.resolution(1),metadata.resolution(2),metadata.resolution(3));
fprintf(fid, 'DimSize = %d %d %d\n',dimensions(1),dimensions(2),dimensions(3));
fprintf(fid, 'ElementType = %s\n',matlabTypeToMHAtype(metadata.datatype));
fprintf(fid, 'ElementDataFile = LOCAL\n');
fwrite(fid,cube,metadata.datatype,'b');

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
       error(['Datatype ' datatype ' not supported by MHA exporter!']);
end
end


