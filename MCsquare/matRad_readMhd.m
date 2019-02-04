function cube = matRad_readMhd(folder,filename)
% matRad mhd file reader
% 
% call
%   cube = matRad_readMhd(folder,filename)
%
% input
%   folder:   folder where the *raw and *mhd file are located
%   filename: filename
%
% output
%   cube:     3D array
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% read header
headerFileHandle = fopen([folder filesep filename],'r');

s = textscan(headerFileHandle, '%s', 'delimiter', '\n');

% read dimensions
idx = find(~cellfun(@isempty,strfind(s{1}, 'DimSize')),1,'first');
dimensions = cell2mat(textscan(s{1}{idx},'DimSize = %f %f %f'));

% read filename of data
idx = find(~cellfun(@isempty,strfind(s{1}, 'ElementDataFile')),1,'first');
tmp = textscan(s{1}{idx},'ElementDataFile = %s');
dataFilename = cell2mat(tmp{1});

% get data type
idx = find(~cellfun(@isempty,strfind(s{1}, 'ElementType')),1,'first');
tmp = textscan(s{1}{idx},'ElementType = MET_%s');
type = lower(cell2mat(tmp{1}));

fclose(headerFileHandle);

%% read data
dataFileHandle = fopen([folder filesep dataFilename],'r');
cube = reshape(fread(dataFileHandle,inf,type),dimensions);
cube = permute(cube,[2 1 3]);
cube = flip(cube,2);
fclose(dataFileHandle);
