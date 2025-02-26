function [cube, metadata] = matRad_readMHD(filename)
% matRad NRRD reader
% 
% call
%   [cube, metadata] = matRad_readMHD(filename)
%
% input
%   filename:   full path to mhd or mha file
%
% output
%   cube:       the read cube
%   metadata:   metadata from header information
%
% References
%   [1] https://itk.org/Wiki/MetaIO/Documentation
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

%% read header
headerFileHandle = fopen(filename,'r');

s = textscan(headerFileHandle, '%s', 'delimiter', '\n');

% read dimensions
idx = find(~cellfun(@isempty,strfind(s{1}, 'DimSize')),1,'first');
dimensions = cell2mat(textscan(s{1}{idx},'DimSize = %f %f %f'));

% Read resolution
idx = find(~cellfun(@isempty,strfind(s{1}, 'ElementSpacing')),1,'first');
tmp = textscan(s{1}{idx},'ElementSpacing = %f %f %f');
resolution = cell2mat(tmp);

%Endian:
idx = find(~cellfun(@isempty,strfind(s{1}, 'BinaryDataByteOrderMSB')),1,'first');
tmp = textscan(s{1}{idx},'BinaryDataByteOrderMSB = %s');
isLittleEndian = cell2mat(tmp{1});
switch isLittleEndian
    case 'True'
        endian = 'b';
    case 'False'
        endian = 'l';
    otherwise 
        matRad_cfg.dispError('Machine format/endian could not be read!');
end


% read filename of data
idx = find(~cellfun(@isempty,strfind(s{1}, 'ElementDataFile')),1,'first');
tmp = textscan(s{1}{idx},'ElementDataFile = %s');
dataFilename = cell2mat(tmp{1});

% Transform Matrix
idx = find(~cellfun(@isempty,strfind(s{1}, 'TransformMatrix')),1,'first');
tmp = textscan(s{1}{idx},'TransformMatrix = %f %f %f %f %f %f %f %f %f');
T = zeros(3);
T(:) = cell2mat(tmp);

if ~isequal(T, diag(ones(1,numel(dimensions))))
    matRad_cfg.dispWarning('Non identity transformation matrix detected in the loaded cube. This might lead to reconstruction inconsistency.')
end
    
% Apply Matlab permutation
% This ensures that the cube is reverted back to the matLab standard
% indexing
Tmatlab = [0 1 0; 1 0 0; 0 0 1];
T = T*Tmatlab;

% get data type
idx = find(~cellfun(@isempty,strfind(s{1}, 'ElementType')),1,'first');
tmp = textscan(s{1}{idx},'ElementType = %s');
type = MHAtypeToMatlab(cell2mat(tmp{1}));

if strcmpi(dataFilename,'LOCAL')
    %read raw block
    test = cast(1,type);
    S = whos('test');
    fseek(headerFileHandle,-S.bytes*prod(dimensions),'eof');
    cube = fread(headerFileHandle,prod(dimensions),type,endian);
    cube = reshape(cube,dimensions);
    cube = permute(cube,abs([1 2 3]*T));
else
    %% read data
    [filepath,~,~] = fileparts(filename);
    dataFileHandle = fopen(fullfile(filepath,dataFilename),'r');
    cube = reshape(fread(dataFileHandle,inf,type,endian),dimensions);
    cube = permute(cube,abs([1 2 3]*T));
    fclose(dataFileHandle);
end
fclose(headerFileHandle);
metadata.resolution = resolution;
metadata.cubeDim = dimensions * T;

end

function newType = MHAtypeToMatlab(datatype)
switch datatype
    case 'MET_FLOAT'
        newType = 'single';     
    case 'MET_DOUBLE'
        newType = 'double';
    case 'MET_UCHAR'
       newType = 'uint8';
    case 'MET_CHAR'
        newType = 'char';
    case 'MET_SHORT'
        newType = 'int16';
    case 'MET_USHORT'
        newType = 'uint16';
    case 'MET_INT'
        newType = 'int32';
    case 'MET_LONG'
        newType = 'int64';
    case 'MET_UINT'
        newType = 'uint32';
    case 'MET_ULONG'
        newType = 'uint64';
    otherwise
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError(['Datatype ' datatype ' not supported by MHD/MHA importer!']);
end
end

