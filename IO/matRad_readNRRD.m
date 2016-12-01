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
nrrdMetaData.fields = cell(0,2);
nrrdMetaData.keys = cell(0,2);
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
            warning(['Could not parse line: "' lineContent '"']);
        elseif isequal(lineContent{1}{2},' ') %space after colon refers to "field"
            nrrdMetaData.fields{end+1,1} = lineContent{1}{1}; %Fieldname
            nrrdMetaData.fields{end,2} = lineContent{1}{3}; %Information
        elseif isequal(lineContent{1}{2},'=') %= after colon refers to key-value pair
            nrrdMetaData.keys{end+1,1} = lineContent{1}{1}; %Key
            nrrdMetaData.keys{end,2} = lineContent{1}{3}; %Value
        else
            warning(['Could not parse line: "' lineContent '"']);
        end        
    end
    currentLine = fgetl(hFile);
end

%% Interpret Headers
%check for the always required data type (and endian if required)
doSwapBytes = false; %If endian differs
typeFieldIx = find(ismember(nrrdMetaData.fields(:,1), 'type'));
if ~isempty(typeFieldIx)
    datatype = getMATLABdataType(nrrdMetaData.fields{typeFieldIx,2});
    %if size of the datatype is > 1 byte, we need endian information
    if ~isequal(datatype(end),'8')
        endianFieldIx = find(ismember(nrrdMetaData.fields(:,1),'endian'));
        if ~isempty(endianFieldIx)
            if ~isequal(nrrdMetaData.fields{endianFieldIx,2},'little') && ~isequal(nrrdMetaData.fields{endianFieldIx,2},'big')
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
            if ~isequal(endian,nrrdMetaData.fields{endianFieldIx,2})
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
dimFieldIx = find(ismember(nrrdMetaData.fields(:,1), 'dimension')); 
if ~isempty(dimFieldIx)
    [metadata.dimension,success] = str2num(nrrdMetaData.fields{dimFieldIx,2});
    if ~success
        error('Could not read required dimension field');
    end
else
    error('Could not find required "dimension" field!');
end

%Check for size / dim length
sizeFieldIx = find(ismember(nrrdMetaData.fields(:,1), 'sizes')); 
if ~isempty(sizeFieldIx)
    sizes = textscan(nrrdMetaData.fields{sizeFieldIx,2},'%d');
    if numel(sizes{1}) ~= metadata.dimension || ~all(sizes{1} > 0) 
        error('Incorrect size definition!');
    end
    metadata.cubeDim = sizes{1}';
else
    error('Could not find required "dimension" field!');
end

%Check for resolution
%Here we need either spacings ore space directions
spacingFieldIx = find(ismember(nrrdMetaData.fields(:,1), 'spacings'));
spaceDirFieldIx = find(ismember(nrrdMetaData.fields(:,1), 'space directions'));
if ~isempty(spacingFieldIx)
    resolutions = textscan(nrrdMetaData.fields{spacingFieldIx,2},'%f');
    if numel(resolutions{1}) ~= metadata.dimension
        error('Incorrect spacings definition');
    end
    metadata.resolution = resolutions{1}';    
    
    %We have default permutation here (required for mapping to MATLAB
    %default order 2 1 3)
    metadata.axisPermutation = [1 2 3];
    
elseif ~isempty(spaceDirFieldIx)
    %space directions are written in vector format
    %first create the scanning string
    vectorstring = '(';
    for c=1:metadata.dimension
        vectorstring = [vectorstring '%f,'];
    end
    vectorstring(end) = ')';
    %Get the vectors
    vectors = textscan(nrrdMetaData.fields{spaceDirFieldIx,2},vectorstring);
    
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
    warning('No Resolution Information available');
end

%find the origin if we have one
originFieldIx = find(ismember(nrrdMetaData.fields(:,1), 'space origin'));
if ~isempty(originFieldIx)
    %first create the scanning string
    vectorstring = '(';
    for c=1:metadata.dimension
        vectorstring = [vectorstring '%f,'];
    end
    originVector = textscan(nrrdMetaData.fields{originFieldIx,2},vectorstring);
    for c=1:metadata.dimension
        metadata.imageOrigin(c) = originVector{c};
    end
end

%Coordinate system - optional
originFieldIx = find(ismember(nrrdMetaData.fields(:,1), 'space'));
if ~isempty(originFieldIx)
    metadata.coordinateSystem = nrrdMetaData.fields{originFieldIx,2};
end
    
%Check for separate file
data_fileFieldIx = find(ismember(nrrdMetaData.fields(:,1), 'data file'));
datafileFieldIx = find(ismember(nrrdMetaData.fields(:,1), 'datafile'));
if ~isempty(data_fileFieldIx) || ~isempty(datafileFieldIx)
    error('Sorry! We currently do not support detached data files!');
    %Proposed workflow:
    %check for data file    
    %close file
    %replace file handle with data file
    %continue without the read part knowing its a new file
end


%% Read Data
%Check for encoding
encodingFieldIx = find(ismember(nrrdMetaData.fields(:,1), 'encoding'));
if isempty(encodingFieldIx)
    error('Could not find required "encoding" field!');
end
switch nrrdMetaData.fields{encodingFieldIx,2}
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
                warning('Java unpacking failed... using temporary files!');
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
    case {'signed char','int8','int8_t'}
        datatype = 'int8';     
    case {'uchar','unsigned char','uint8','uint8_t'}
        datatype = 'uint8';        
    case {'short','short int','signed short','signed short int','int16','int16_t'}
        datatype = 'int16';        
    case {'ushort','unsigned short','unsigned short int','uint16','uint16_t'}
        datatype = 'uint16';        
    case {'int','signed int','int32','int32_t'}
        datatype = 'int32';        
    case {'uint','unsigned int','uint32','uint32_t'}
        datatype = 'uint32';        
    case {'longlong','long long','long long int','signed long long','signed long long int','int64','int64_t'}
        datatype = 'int64';        
    case {'ulonglong', 'unsigned long long', 'unsigned long long int','uint64','uint64_t'}
        datatype = 'uint64';        
    case 'float'
        datatype = 'single';        
    case 'double'
        datatype = 'double';        
    otherwise
        error('Could not identify datatype for NRRD data');
end

end

