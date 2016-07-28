function [cube, metadata] = matRad_readNRRD(filename)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad NRRD reader
% 
% call
%   matRad_readNRRD(filename)
%
% input
%   filename:   full path to nrrd file
%
% output
%   cube:       the read cube
%   metadata:   metadata from header information
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

% Open file.
hFile = fopen(filename, 'r');
if hFile <= 0
    error('Could not open NRRD file!');
end
cleaner = onCleanup(@() fclose(hFile));

%% Determine NRRD Version
nrrdLine = fgetl(hFile);
regTokens = regexp(nrrdLine,'NRRD00(?:\.)?0([1-5])','tokens');
if isempty(regTokens)
    error('Invalid Header line! Could not identify NRRD version!');
end
nrrdVersion = str2num(regTokens{1}{1});

if nrrdVersion > 5
    error('NRRD version > 5 not supported!');
end

%% Read header
nrrdMetaData.fields = struct();
nrrdMetaData.keys = struct();
nrrdMetaData.comments = cell(0);

currentLine = fgetl(hFile);
while ~isempty(currentLine) && ischar(currentLine) %NRRD separates data from header by an empty line
    %Check for comment
    if isequal(currentLine(1),'#')
        nrrdMetaData.comments{end+1} = currentLine;
    else
        %Parse the line
        lineContent = regexp(currentLine, '(.+):(=|\s)(.+)', 'tokens');
        if isempty(lineContent)
            warn(['Could not parse line: "' lineContent '"']);
        elseif isequal(lineContent{1}{2},' ') %space after colon refers to "field"
            fieldname = regexprep(lineContent{1}{1},'[-/\s]','_');
            nrrdMetaData.fields.(fieldname) = lineContent{1}{3};
        elseif isequal(lineContent{1}{2},'=') %= after colon refers to key-value pair
            key = regexprep(lineContent{1}{1},'[-/\s]','_');
            nrrdMetaData.keys.(key) = lineContent{1}{3};
        else
            warn(['Could not parse line: "' lineContent '"']);
        end        
    end
    currentLine = fgetl(hFile);
end

%% Interpret Headers
%check for the always required data type (and endian if required)
doSwapBytes = false; %If endian differs
if isfield(nrrdMetaData.fields,'type')
    datatype = getMATLABdataType(nrrdMetaData.fields.type);
    %if size of the datatype is > 1 byte, we need endian information
    if ~isequal(datatype(end),'8')
        if isfield(nrrdMetaData.fields,'endian')
            if ~isequal(nrrdMetaData.fields.endian,'little') && ~isequal(nrrdMetaData.fields.endian,'big')
                error(['Datatype is ' datatype ', thus endian information is required but could not be interpreted!']);
            end;
            %Now we compare the file endian to the system endian
            %First acquire system endian
            [~,~,endian] = computer();
            if isequal(endian,'B')
                endian = 'big';
            end;
            if isequal(endian,'L')
                endian = 'little';
            end;
            %now compare to file endian and set flag if appropriate
            if ~isequal(endian,nrrdMetaData.fields.endian)
                doSwapBytes = true;
            end
        else
            error(['Datatype is ' datatype ', thus endian information is required but could not be found!']);
        end
    end
    metadata.datatype = datatype;
    
else
    error('Could not find required "type" field!');
end

%Check for the always required image dimension
if isfield(nrrdMetaData.fields,'dimension')
    [metadata.dimension,success] = str2num(nrrdMetaData.fields.dimension);
    if ~success
        error('Could not read required dimension field');
    end
else
    error('Could not find required "dimension" field!');
end

%Check for size / dim length
if isfield(nrrdMetaData.fields,'sizes')
    sizes = textscan(nrrdMetaData.fields.sizes,'%d');
    if numel(sizes{1}) ~= metadata.dimension || ~all(sizes{1} > 0) 
        error('Incorrect size definition!');
    end
    metadata.cubeDim = sizes{1}';
else
    error('Could not find required "dimension" field!');
end

%Check for resolution
%Here we need either spacings ore space directions
if isfield(nrrdMetaData.fields,'spacings')
    resolutions = textscan(nrrdMetaData.fields.spacings,'%f');
    if numel(resolutions{1}) ~= metadata.dimension
        error('Incorrect spacings definition');
    end
    metadata.resolution = resolutions{1}';    
    
    %We have default permutation here (required for mapping to MATLAB
    %default order 2 1 3)
    metadata.axisPermutation = [1 2 3];
    
elseif isfield(nrrdMetaData.fields,'space_directions')
    %space directions are written in vector format
    %first create the scanning string
    vectorstring = '(';
    for c=1:metadata.dimension
        vectorstring = [vectorstring '%f,'];
    end
    vectorstring(end) = ')';
    %Get the vectors
    vectors = textscan(nrrdMetaData.fields.space_directions,vectorstring);
    
    %At the moment we only support cartesian basis vectors
    %this gives us the permutation to align with the MATLAB default
    %ordering of 2 1 3
    for c=1:metadata.dimension
        %check if cartesian basis vector        
        currentAxis = find(vectors{c});
        
        if numel(find(vectors{c})) ~= 1
            error('Sorry! We currently only support spaces with cartesian basis!'); 
        end        
        metadata.axisPermutation(c) = currentAxis*sign(vectors{c}(currentAxis));
        metadata.resolution(c) = vectors{c}(currentAxis);       
    end
   
else
    warn('No Resolution Information available');
end

if isfield(nrrdMetaData.fields,'space_origin')
    %first create the scanning string
    vectorstring = '(';
    for c=1:metadata.dimension
        vectorstring = [vectorstring '%f,'];
    end
    originVector = textscan(nrrdMetaData.fields.space_origin,vectorstring);
    for c=1:metadata.dimension
        metadata.imageOrigin(c) = originVector{c};
    end
end

%Coordinate system - optional
if isfield(nrrdMetaData.fields,'space')
    metadata.coordinateSystem = nrrdMetaData.fields.space;
end
    
%Check for separate file
if isfield(nrrdMetaData.fields,'data file') || isfield(nrrdMetaData.fields,'datafile')
    error('Sorry! We currently do not support detached data files!');
    %Proposed workflow:
    %check for data file    
    %close file
    %replace file handle with data file
    %continue without the read part knowing its a new file
end


%% Read Data
%Check for encoding
if ~isfield(nrrdMetaData.fields,'encoding')
    error('Could not find required "encoding" field!');
end
switch nrrdMetaData.fields.encoding
    case 'raw'
        cube = fread(hFile,prod(metadata.cubeDim), metadata.datatype);
    case {'txt','text','ascii'}
        cube = cast(fscanf(hFile,'%f'),metadata.datatype);
    case 'hex'
        error('Sorry: NRRD hex file not yet supported!');
    case {'gz','gzip'}
        compressedByteArray = fread(hFile, inf, 'uint8');
        
        javaUnpackSuccessful = false;
        
        if isempty(javachk('jvm'))
            %Unzip with java
            try
                javaByteStream = java.io.ByteArrayInputStream(compressedByteArray);
                gzipStream = java.util.zip.GZIPInputStream(javaByteStream);
                javaByteOutputStream = java.io.ByteArrayOutputStream();
                org.apache.commons.io.IOUtils.copy(gzipStream,javaByteOutputStream);
                gzipStream.close();
                cube = typecast(javaByteOutputStream.toByteArray(),metadata.datatype)';
                javaByteStream.close();
                javaByteOutputStream.close();
                javaUnpackSuccessful = true;
            catch
                warn('Java unpacking failed... using temporary files!');
            end
        end
        if ~javaUnpackSuccessful
            %Copy content to temporary file
            tmpName = tempname();
            tmpFile = [tmpName '.gz'];
            hFileTmp = fopen(tmpFile, 'wb');
            
            if hFileTmp <= 0
                error('Could not open temporary file for GZIP!');
            end
            
            fwrite(hFileTmp, compressedByteArray, 'uint8');
            fclose(hFileTmp);
            
            %Unzip with gzip
            gunzip(tmpFile);
            
            %Read the uncompressed file
            hFileTmp = fopen(tmpName, 'rb');
            if hFileTmp <= 0
                error('Could not open unpacked file!');
            end
            cleanTmpFile = onCleanup(@() fclose(hFileTmp));

            cube = fread(hFileTmp, prod(metadata.cubeDim), metadata.datatype);
        end
    case {'bz2','bzip2'}
        error('Sorry: bzip compression not yet supported!');
    otherwise 
        error(['Undefined NRRD encoding scheme: ' nrrdMetaData.encoding]);
end

%maybe we need to correct the byte ordering (endian)
if doSwapBytes
    swapbytes(cube);
end

%first we shape the data into a cube
cube = reshape(cube,metadata.cubeDim);

%now we have to do the permutations, 2 1 3 ... is the MATLAB default
%We create a transform matrix that transforms a permutation to MATLAB
%default
permutationTransformMatrix = diag(ones(metadata.dimension,1));
permutationTransformMatrix(1:2,1:2) = flip(permutationTransformMatrix(1:2,1:2));

applyPermutation = permutationTransformMatrix*metadata.axisPermutation';

cube = permute(cube,applyPermutation);

end

%This function wraps datatype definitions to the one used by MATLAB
function datatype = getMATLABdataType(typestring)
%The typedefinitions are taken directly from Section 5 of the NRRD
%definition
switch typestring
    case {'signed char', 'int8', 'int8_t'}
        datatype = 'int8';     
    case {'uchar', 'unsigned char', 'uint8', 'uint8_t'}
        datatype = 'uint8';        
    case {'short', 'short int', 'signed short', 'signed short int', ...
            'int16', 'int16_t'}
        datatype = 'int16';        
    case {'ushort', 'unsigned short', 'unsigned short int', 'uint16', ...
            'uint16_t'}
        datatype = 'uint16';        
    case {'int', 'signed int', 'int32', 'int32_t'}
        datatype = 'int32';        
    case {'uint', 'unsigned int', 'uint32', 'uint32_t'}
        datatype = 'uint32';        
    case {'longlong', 'long long', 'long long int', 'signed long long', ...
            'signed long long int', 'int64', 'int64_t'}
        datatype = 'int64';        
    case {'ulonglong', 'unsigned long long', 'unsigned long long int', ...
            'uint64', 'uint64_t'}
        datatype = 'uint64';        
    case 'float'
        datatype = 'single';        
    case 'double'
        datatype = 'double';        
    otherwise
        error('Could not identify datatype for NRRD data');
end

end

