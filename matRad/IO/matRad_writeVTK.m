function matRad_writeVTK(filepath,cube,metadata)
% matRad function to write vtk cubes
% 
% call
%   matRad_writeVTK(filepath,cube,metadata)
%
% input
%   filepath:   full filename (with extension)
%   cube:       3D array to be written into file
%   metadata:   struct of metadata. Writer will wrap the existing metadata 
%               to VTK standard-specific fields 
%               Necessary fieldnames are:
%               - resolution: [x y z]
%               - datatype: numeric MATLAB-Datatype
%
% output
%   file will be written to disk
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

% Sanity checks and restrictions
dimensions = size(cube);
if numel(dimensions) ~= 3
    error('Sorry! matRad only supports 3-dimensional VTK output');
end

fid = fopen(filepath, 'wb');
if fid <= 0
    error('Could not open VTK destination file!');
end
cleaner = onCleanup(@() fclose(fid));

%We perform 
if isfield(metadata,'axisPermutation')
    cube = permute(cube,metadata.axisPermutation);
end


fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'vtk output\n');
fprintf(fid, 'BINARY\n');
fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', dimensions(1),dimensions(2),dimensions(3));
fprintf(fid, 'SPACING    %f   %f   %f\n',metadata.resolution(1),metadata.resolution(2),metadata.resolution(3));
fprintf(fid, 'ORIGIN    %f   %f   %f\n',metadata.imageOrigin(1),metadata.imageOrigin(2),metadata.imageOrigin(3));
fprintf(fid, 'POINT_DATA   %d\n',prod(dimensions));
if isfield(metadata,'dataName')
    dataName = metadata.dataName;
else
    dataName = 'scalars';
end
fprintf(fid, 'SCALARS %s %s\n',dataName,metadata.datatype);
fprintf(fid, 'LOOKUP_TABLE default\n');
fwrite(fid,cube,metadata.datatype,'b');

end


