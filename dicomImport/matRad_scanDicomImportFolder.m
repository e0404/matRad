function [ fileList, patientList ] = matRad_scanDicomImportFolder( patDir )
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to scan a folder for dicom data
% 
% call
%   [ fileList, patientList ] = matRad_scanDicomImportFolder( patDir )
%
% input
%   patDir:         folder to be scanned
%
% output
%   fileList:       matlab struct with a list of dicom files including meta
%                   infomation (type, series number etc.)
%   patientList:    list of patients with dicom data in the folder
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

%% get all files in search directory

% get information about main directory
mainDirInfo = dir(patDir);
% get index of subfolders
dirIndex = [mainDirInfo.isdir];
% list of filenames in main directory
fileList = {mainDirInfo(~dirIndex).name}';
patientList = 0;

% create full path for all files in main directory
if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(patDir,x),...
        fileList, 'UniformOutput', false);
    
    %% check for dicom files and differentiate patients, types, and series
    numOfFiles = numel(fileList(:,1));
    h = waitbar(0,'Please wait...');
    %h.WindowStyle = 'Modal';
    steps = numOfFiles;
    for i = numOfFiles:-1:1
        waitbar((numOfFiles+1-i) / steps)
        try % try to get DicomInfo
            info = dicominfo(fileList{i});
        catch
            fileList(i,:) = [];
            continue;
        end
        try
            fileList{i,2} = info.Modality;
        catch
            fileList{i,2} = NaN;
        end
        try
            fileList{i,3} = info.PatientID;
        catch
            fileList{i,3} = NaN;
        end
        try
            fileList{i,4} = info.SeriesInstanceUID;
        catch
            fileList{i,4} = NaN;
        end
        try
            fileList{i,5} = num2str(info.SeriesNumber);
        catch
            fileList{i,5} = NaN;
        end
        try
            fileList{i,6} = info.PatientName.FamilyName;
        catch
            fileList{i,6} = NaN;
        end
        try
            fileList{i,7} = info.PatientName.GivenName;
        catch
            fileList{i,7} = NaN;
        end
        try
            fileList{i,8} = info.PatientBirthDate;
        catch
            fileList{i,8} = NaN;
        end
        try
            if strcmp(info.Modality,'CT')
                fileList{i,9} = num2str(info.PixelSpacing(1));
            else
                fileList{i,9} = NaN;
            end
        catch
            fileList{i,9} = NaN;
        end
        try
            if strcmp(info.Modality,'CT')
                fileList{i,10} = num2str(info.PixelSpacing(2));
            else
                fileList{i,10} = NaN;
            end
        catch
            fileList{i,10} = NaN;
        end
        try
            if strcmp(info.Modality,'CT')
                fileList{i,11} = num2str(info.SliceThickness);
            else
                fileList{i,11} = NaN;
            end
        catch
            fileList{i,11} = NaN;
        end
        
        matRad_progress(numOfFiles+1-i, numOfFiles);
        
    end
    close(h)
    
    if ~isempty(fileList)
        patientList = unique(fileList(:,3));
        % check if there is at least one RT struct and one ct file
        % available per patient
        for i = numel(patientList):-1:1
            if sum(strcmp('CT',fileList(:,2)) & strcmp(patientList{i},fileList(:,3))) < 1 || ...
               sum(strcmp('RTSTRUCT',fileList(:,2)) & strcmp(patientList{i},fileList(:,3))) < 1
                patientList(i) = [];
            end
        end
        
        if isempty(patientList)
            msgbox('No patient found with DICOM CT _and_ RT structure set in patient directory!', 'Error','error');
        end
    else
        msgbox('No DICOM files found in patient directory!', 'Error','error');
        %h.WindowStyle = 'Modal';
        %error('No DICOM files found in patient directory');
    end
else
    h = msgbox('Search folder empty!', 'Error','error');
    %h.WindowStyle = 'Modal';
    %error('Search folder empty')
    
end


end

