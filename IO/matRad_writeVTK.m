function matRad_writeVTK(filepath,cube,metadata)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to write vtk cubes
% 
% call
%   matRad_writeVTK(cube,metadata,filename)
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sanity checks and restrictions
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


