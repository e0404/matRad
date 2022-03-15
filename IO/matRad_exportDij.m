function matRad_exportDij(filename,dij,stf,metadata)
% matRad physical dose writer
%
% call
%   matRad_exportDij(filename,dij,stf,...
%                    additionalFields,additionalKeyValuePairs)
%
% input
%   filename:   full output path, including the file extension
%   dij:        matRad dij struct
%   stf:        matRad stf struct
%   metadata:   struct of metadata
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

matRad_cfg = MatRad_Config.instance();

if nargin<4
    metadata = struct();
end


%% Prepare Metadata

if ~isfield(metadata,'delimiter')
    metadata.delimiter = '\t'; %Default delimiter
end

if ~isfield(metadata,'numScen')
    metadata.numScen = 1; %Default scenario
end

if ~isfield(metadata,'individualFiles')
    metadata.individualFiles = false; %Default individual files option
end

if ~isfield(metadata,'extension')
    
    lastdot_pos = find(filename == '.', 1, 'last');
    extension = filename(lastdot_pos+1:end);
    
    if strcmp(extension,'txt') || strcmp(extension,'bin')
        metadata.extension = extension; %Default fileType
    else
        metadata.extension = 'txt'; %Default fileType
    end
    
end

%% Setup Header

header = sprintf('# %s %s\n',metadata.extension,'file');

%add matRad specific comment
header = header_addComment(header,'Created With matRad - An open source multi-modality radiation treatment planning sytem');

%% Write File
try
    
    %Set up parent export folder and full file path
    if ~(isfolder('dijExport'))
        mkdir(matRad_cfg.matRadRoot, 'dijExport');
    end
    
    folderPath = [matRad_cfg.matRadRoot filesep 'dijExport' filesep];
    
    if metadata.individualFiles
        
        totalNumOfBixels = 1;
        
        %Create a file for each beam
        for i = 1:dij.numOfBeams
            
            %Set a filename for i-th beam file
            lastdot_pos = find(filename == '.', 1, 'last');
            
            filename_ith = filename(1:lastdot_pos-1);
            filename_ith = [filename_ith '_' num2str(i)];
            
            %Add gantryAngle field to i-th beam header
            header_ith = header_addIntField(header,'gantry angle',stf(i).gantryAngle);
            %Add couchAngle field to i-th beam header
            header_ith = header_addIntField(header_ith,'couch angle', stf(i).couchAngle);
            %Add totalNumOfBixels field to i-th beam header
            header_ith = header_addIntField(header_ith,'total number of bixels', stf(i).totalNumOfBixels);
            %Add dimensions of dose grid field to i-th beam header
            dimensions = strcat(num2str(dij.doseGrid.dimensions(1)),'|',num2str(dij.doseGrid.dimensions(2)),'|',num2str(dij.doseGrid.dimensions(3)));
            header_ith = header_addStringField(header_ith,'dose grid dimensions', dimensions);
            %Add column headers
            header_ith = header_addComment(header_ith,'voxelID bixelID physicalDose[Gy]');
            
            %
            numOfBixels = stf(i).totalNumOfBixels;
            
            %Read physical dose from non zeros dij
            [ix,iy,vals] = find(dij.physicalDose{metadata.numScen}(:,totalNumOfBixels:totalNumOfBixels+numOfBixels-1));
            data=zeros(nnz(vals),3);
            data(:,1) = ix;
            data(:,2) = iy+totalNumOfBixels-1;
            data(:,3) = vals;
            
            if strcmp(metadata.extension,'txt')
                
                %Write Header to file with the separating blank line to i-th beam
                fileHandle = fopen([folderPath filename_ith '.' metadata.extension],'w');
                fprintf(fileHandle,'%s\n',header_ith);
                
                %Append data to file to i-th beam
                %writematrix(data,filename_tmp,'Delimiter',metadata.delimiter,'-append'); % If you use r2019b matlab version
                dlmwrite([filename_ith '.' metadata.extension],data,'delimiter',metadata.delimiter,'-append');
                
                fclose(fileHandle);
                
            elseif strcmp(metadata.extension,'bin')
                
                %Append data to file to i-th beam
                fileHandle = fopen([folderPath filename_ith '.' metadata.extension],'w');
                fwrite(fileHandle,uint32(ix),'uint32');
                fwrite(fileHandle,uint32(iy),'uint32');
                fwrite(fileHandle,vals,'double');
                fclose(fileHandle);
                
                %Write an additional header file
                headerHandle = fopen([folderPath filename_ith '_header.txt'],'w');
                fprintf(headerHandle,'%s\n',header_ith);
                fclose(headerHandle);
                
            end
            
            totalNumOfBixels=totalNumOfBixels+numOfBixels;
            
        end
        
    else
        
        %Add info about each beam
        for i = 1:dij.numOfBeams
            %Add info about i-th beam
            header = header_addIntField(header,'Beam',i);
            %Add gantryAngle field to i-th beam header
            header = header_addIntField(header,'gantry angle',stf(i).gantryAngle);
            %Add couchAngle field to i-th beam header
            header = header_addIntField(header,'couch angle', stf(i).couchAngle);
            %Add totalNumOfBixels field to i-th beam header
            header = header_addIntField(header,'total number of bixels', stf(i).totalNumOfBixels);
            %Add dimensions of dose grid field to header
            dimensions = strcat(num2str(dij.doseGrid.dimensions(1)),'|',num2str(dij.doseGrid.dimensions(2)),'|',num2str(dij.doseGrid.dimensions(3)));
            header = header_addStringField(header,'dose grid dimensions', dimensions);
        end
        
        %Set a filename
        filename = filename(1:lastdot_pos-1);
        
        %Add column headers
        header = header_addComment(header,'voxelID bixelID physicalDose[Gy]');
        
        %Read physical dose from non zeros dij
        [ix,iy,vals] = find(dij.physicalDose{metadata.numScen});
        data(:,1) = ix;
        data(:,2) = iy;
        data(:,3) = vals;
        
        if strcmp(metadata.extension,'txt')
            
            %Write Header to file with the separating blank line to i-th beam
            fileHandle = fopen([folderPath filename '.' metadata.extension],'w');
            fprintf(fileHandle,'%s\n',header);
            
            %Append data to file
            %writematrix(data,filename,'Delimiter',metadata.delimiter,'-append'); % If you use r2019b matlab version
            dlmwrite([folderPath filename '.' metadata.extension],data,'delimiter',metadata.delimiter,'-append');
            
            fclose(fileHandle);
            
        elseif strcmp(metadata.extension,'bin')
            
            %Append data to file
            fileHandle = fopen([folderPath filename '.' metadata.extension],'w');
            fwrite(fileHandle,uint32(ix),'uint32');
            fwrite(fileHandle,uint32(iy),'uint32');
            fwrite(fileHandle,vals,'double');
            fclose(fileHandle);
            
            %Write an additional header file
            headerHandle = fopen([folderPath filename '_header.txt'],'w');
            fprintf(headerHandle,'%s\n',header);
            fclose(headerHandle);
            
        end
        
    end
    
catch MExc
    %if something failed while writing, close all files and display error
    fclose('all');
    fprintf(2,'File %s could not be written!\n',filename);
    if(matRad_cfg.isOctave)
        error(MExc);
    else
        throw(MExc);  
    end
end

fprintf(1,'Dij exported successfully into %s.\n',strcat(folderPath,filename,'.',metadata.extension));

%Used to add comments to the header
    function newHeader = header_addComment(header,comment)
        newHeader = sprintf('%s# %s\n',header,comment);
    end

%Used to add int fields to the header
    function newHeader = header_addIntField(header,fieldName,fieldValue)
        newHeader = sprintf('%s# %s: %d\n',header,fieldName,fieldValue);
    end

%Used to add string fields to the header
    function newHeader = header_addStringField(header,fieldName,fieldValue)
        newHeader = sprintf('%s# %s: %s\n',header,fieldName,fieldValue);
    end

end
