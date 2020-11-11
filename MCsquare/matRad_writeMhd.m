function matRad_writeMhd(cube,resolution,filename)


% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% write header file
fileHandle = fopen(filename,'w');

fprintf(fileHandle,'ObjectType = Image\n');
fprintf(fileHandle,'NDims = 3\n');
fprintf(fileHandle,'BinaryData = True\n');
fprintf(fileHandle,'BinaryDataByteOrderMSB = False\n');
fprintf(fileHandle,'CompressedData = False\n');
fprintf(fileHandle,'TransformMatrix = 1 0 0 0 1 0 0 0 1\n');
fprintf(fileHandle,'Offset = 0 0 0\n');
fprintf(fileHandle,'CenterOfRotation = 0 0 0\n');
fprintf(fileHandle,'AnatomicalOrientation = RAI\n');
fprintf(fileHandle,'ElementSpacing = %f %f %f\n',resolution);
fprintf(fileHandle,'DimSize = %d %d %d\n',size(cube,2),size(cube,1),size(cube,3));
fprintf(fileHandle,'ElementType = MET_DOUBLE\n');
filenameRaw = [filename(1:end-4) '.raw'];
fprintf(fileHandle,'ElementDataFile = %s\n',filenameRaw);

fclose(fileHandle);

%% write data file
dataFileHandle = fopen(filenameRaw,'w');

cube = flip(cube,2);
cube = permute(cube,[2 1 3]);

fwrite(dataFileHandle,cube(:),'double');
fclose(dataFileHandle);
