function matRad_writeNRRD(filename,cube,metadata)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad NRRD writer
% 
% call
%   matRad_writeNRRD(filename,cube,datatype,...
%                    additionalFields,additionalKeyValuePairs)
%
% input
%   filename:   full output path, including the nrrd extension
%   cube:       cube that is to be written
%   metadata:   struct of metadata. Writer will wrap the existing metadata 
%               to nrrd standard-specific fields [1].
%
% References
%   [1] http://teem.sourceforge.net/nrrd/format.html5
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
%% Setup Header
%The NRRD header description line
version = 5;
nrrdVersionStringPrefix = 'NRRD000';
nrrdVersionString = [nrrdVersionStringPrefix num2str(version)];

header = sprintf('%s\n',nrrdVersionString);

%add matRad specific comment
header = header_addComment(header,'Created With matRad - An open source multi-modality radiation treatment planning sytem');

%Add Datatype field
header = header_addField(header,'type',matlabTypeToNRRD(metadata.datatype));
%NRRD wants to know how the bytes are aligned for datatypes > 1 byte
%(little endian vs big endian)
try
    z = zeros(1,metadata.datatype);
    varInfo = whos('z');
    if varInfo.bytes > 1
        [~,~,endian] = computer;
        switch endian
            case 'L'
                header = header_addField(header,'endian','little');
            case 'B'
                header = header_addField(header,'endian','big');
            otherwise
                error('Unknown endian!');
        end
    end
catch
    error(['Unknown datatype: ' metadata.datatype']);
end


%Check if compression was defined
if ~isfield(metadata,'compress')
    metadata.compress = false;
end  
%If we use compression, try to compress, otherwise use raw data 
if metadata.compress
    try
        %We write a temporary file to gzip with matlab
        %It might be also possible to do this with java directly as stream
        tmpRawFile = [tempname() '.bin'];
        fTmp = fopen(tmpRawFile, 'wb');
        
        fwrite(fTmp,cube,metadata.datatype);
        fclose(fTmp);
        
        %Now zip the file
        fNameZip = gzip(tmpRawFile);
        fNameZip = fNameZip{1};
        
        %We need this to determine the file size
        fDir = dir(fNameZip);
        
        %Read the zipped data
        fTmp = fopen(fNameZip, 'rb');
        fileData = fread(fTmp,fDir.bytes,'uint8');
        fclose(fTmp);
        %Delete temporary files
        delete(tmpRawFile,fNameZip);
        
        %Set the datatype to binary for writing the zipped data
        metadata.datatype = 'uint8';
        header = header_addField(header,'encoding','gzip');
    catch
        warn('Could not open temporary file, writing without compression!');
        header = header_addField(header,'encoding','raw');
        fileData = cube;
    end
else
    header = header_addField(header,'encoding','raw');
    fileData = cube;
end


%Dimensionality
cubeDim = size(cube);
header = header_addField(header,'dimension',num2str(numel(cubeDim)));
strCubeDim = vec2str(cubeDim);
header = header_addField(header,'sizes',strCubeDim);

%Coordinate system - optional
if isfield(metadata,'coordinateSystem')
    header = header_addField(header,'space',metadata.coordinateSystem);
end

%Reference Point - optional
if isfield(metadata,'imageOrigin')
    vecOrigin = matVec2nrrdVec(metadata.imageOrigin);
    header = header_addField(header,'space origin',vecOrigin);
end

%Axis permutation 
if isfield(metadata,'resolution')
    %If we have more space information, write the space axis
    if isfield(metadata,'axisPermutation')
        strDirection = '';
        for dim = 1:numel(cubeDim)
            dimVec = zeros(1,numel(cubeDim));
            dimVec(dim) = 1;
            dimVec = dimVec(metadata.axisPermutation);
            dimVec = dimVec .* metadata.resolution;
            strDirection = [strDirection matVec2nrrdVec(dimVec) ' '];
        end
        header = header_addField(header,'space directions',strDirection);
    %Otherwise only write resolution
    else
        strSpacings = vec2str(metadata.resolution);
        header = header_addField(header,'spacings',strSpacings);
    end
end       

%Additional key-value-pairs if present
%Will be added in a later release
%{
if nargin > 3
    keys = fieldnames(additionalKeyValuePairs);
    for f = 1:numel(fields)
        key = keys{f};
        value = additionalKeyValuePairs.(keys{f});
        
        %Format the description to fit
        if isvector(description)
            value = mat2str(description);
            value = description(2:end-1);
        elseif isnumeric(description)
            value = num2str(description);
        else 
     
        end
        
        if and(~ischar(value),~iscellstr(value))
            warning(['Cannot convert value for key ' key ' to string, key-value-pair will be ignored!']);
            continue;
        end
        
        header = header_addKeyValuePair(header,key,value);
    end;
end
%}

%% Write File
try
    %Write Header to file with the separating blank line
    fileHandle = fopen(filename,'w');
    fprintf(fileHandle,'%s\n',header);
    
    %Write raw or compressed data to file
    fwrite(fileHandle,fileData,metadata.datatype);
    fclose(fileHandle);
    
catch MExc
    fclose('all');
    error(sprintf('File %s could not be written!\n%s',filename,getReport(MExc)));
end

end

%Used to add comments to the header
function newHeader = header_addComment(header,comment)
    newHeader = sprintf('%s# %s\n',header,comment);
end

%Used to add fields to the header
function newHeader = header_addField(header,fieldName,fieldDescription)
    newHeader = sprintf('%s%s: %s\n',header,fieldName,fieldDescription);
end

%Used to add key-value pairs to the header
function newHeader = header_addKeyValuePair(header,key,value)
    value = strrep(value,':=','');
    newHeader = sprintf('%s%s:=%s\n',header,key,value);
end

%Used to format dimenional info according to nrrd specification
function strOutput = vec2str(V)
    strOutput = mat2str(V);
    strOutput = strOutput(2:end-1);
end

%Used to format a matlab vector to nrrd vector with brackets
function strOutput = matVec2nrrdVec(V)
    strOutput = '(';
    for v = 1:numel(V)
        strOutput = [strOutput num2str(V(v)) ','];
    end
    strOutput = [strOutput(1:end-1) ')'];
end

%Used to map matlab definitions of datatypes to the NRRD specs
%seems to be only single has to be converted to float
function newType = matlabTypeToNRRD(datatype)
switch datatype
    case 'single'
        newType = 'float';
    otherwise
        newType = datatype;
end
end
        
