function matRad_writeNRRD(filename,cube,datatype,additionalFields,additionalKeyValuePairs)
%MATRAD_WRITENRRD Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    error('Filename and cube must be specified');
end

if nargin < 3
    datatype = 'double'; %default datatype
    fprintf('No datatype specified, using %s\n',datatype);
end
    
%% Setup Header
%The NRRD header description line
version = 3;
nrrdVersionStringPrefix = 'NRRD000';
nrrdVersionString = [nrrdVersionStringPrefix num2str(version)];

header = sprintf('%s\n',nrrdVersionString);

%add matRad specific comment
header = header_addComment(header,'Created With matRad - An open source multi-modality radiation treatment planning sytem');

%Add Datatype field
header = header_addField(header,'type',datatype);
header = header_addField(header,'encoding','raw');
if or(strcmp(datatype,'float'),strcmp(datatype,'double'))
    header = header_addField(header,'endian','little');
end
%Dimensionality
cubeDim = size(cube);
header = header_addField(header,'dimension',num2str(numel(cubeDim)));
cubeDimString = mat2str(cubeDim);
cubeDimString = cubeDimString(2:end-1); %Remove Brackets
header = header_addField(header,'sizes',cubeDimString);

%Additional fields if present
if nargin >= 3
    fields = fieldnames(additionalFields);
    for f = 1:numel(fields)
        field = fields{f};
        description = additionalFields.(fields{f});
        
        %Format the description to fit
        if isvector(description)
            description = mat2str(description);
            description = description(2:end-1);
        elseif isnumeric(description)
            description = num2str(description);
        else 
            
        end
        
        header = header_addField(header,field,description);
    end;
end    

%% Prepare Data
%Permute Dimensions since matRad swaps the x and y dims
cube = permute(cube,[2 1 3]);

%% Write File
try
    %Write Header to file with the separating blank line
    fileHandle = fopen(filename,'w');
    fprintf(fileHandle,'%s\n',header);

    %Write data to file
    fwrite(fileHandle,cube,datatype);

    fclose(fileHandle);
catch MExc
    fclose('all');
    error(sprintf('File %s could not be written!\n%s',filename,getReport(MExc)));
end
fprintf('File written to %s...\n',filename);

end

function newHeader = header_addComment(header,comment)
    newHeader = sprintf('%s# %s\n',header,comment);
end

function newHeader = header_addField(header,fieldName,fieldDescription)
    newHeader = sprintf('%s%s: %s\n',header,fieldName,fieldDescription);
end

function newHeader = header_addKeyValuePair(header,key,value)
    newHeader = sprintf('%s%s:=%s\n',header,key,value);
end

