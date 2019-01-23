function cube = matRad_readMhd(folder,filename)


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
